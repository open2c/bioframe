from __future__ import division, print_function
from collections import OrderedDict
import tempfile
import json
import six
import os

import numpy as np
import pandas as pd
import pyfaidx
import pysam
from ..core import argnatsort
from ..schemas import SCHEMAS


def read_table(filepath_or, schema=None, **kwargs):
    kwargs.setdefault('sep', '\t')
    kwargs.setdefault('header', None)
    if isinstance(filepath_or, six.string_types) and filepath_or.endswith('.gz'):
        kwargs.setdefault('compression', 'gzip')
    if schema is not None:
        try:
            kwargs.setdefault('names', SCHEMAS[schema])
        except (KeyError, TypeError):
            if isinstance(schema, six.string_types):
                raise ValueError("TSV schema not found: '{}'".format(schema))
            kwargs.setdefault('names', schema)
    return pd.read_csv(filepath_or, **kwargs)


def read_chromsizes(filepath_or,
                    name_patterns=(r'^chr[0-9]+$', r'^chr[XY]$', r'^chrM$'),
                    natsort=True,
                    **kwargs):
    """
    Parse a ``<db>.chrom.sizes`` or ``<db>.chromInfo.txt`` file from the UCSC
    database, where ``db`` is a genome assembly name.

    Parameters
    ----------
    filepath_or : str or file-like
        Path or url to text file, or buffer.
    name_patterns : sequence, optional
        Sequence of regular expressions to capture desired sequence names.
        Each corresponding set of records will be sorted in natural order.
    all_names : bool, optional
        Whether to return all contigs listed in the file. Default is
        ``False``.

    Returns
    -------
    Series of integer bp lengths indexed by sequence name.

    See also
    --------
    * UCSC assembly terminology: <http://genome.ucsc.edu/FAQ/FAQdownloads.html#download9>
    * NCBI assembly terminology: <https://www.ncbi.nlm.nih.gov/grc/help/definitions>

    """
    if kwargs.pop('all_names', False):
        name_patterns = 'all'
    if isinstance(filepath_or, six.string_types) and filepath_or.endswith('.gz'):
        kwargs.setdefault('compression', 'gzip')

    chromtable = pd.read_csv(
        filepath_or, sep='\t', usecols=[0, 1],
        names=['name', 'length'], dtype={'name':str}, **kwargs)
    
    if name_patterns != 'all':
        parts = []
        for pattern in name_patterns:
            if not len(pattern): continue
            part = chromtable[chromtable['name'].str.contains(pattern)]
            if natsort:
                part = part.iloc[argnatsort(part['name'])]
            parts.append(part)
        chromtable = pd.concat(parts, axis=0)
    
    chromtable.index = chromtable['name'].values
    return chromtable['length']


def read_gapfile(filepath_or_fp, chroms=None, **kwargs):
    gap = pd.read_csv(
        filepath_or_fp,
        sep='\t',
        names=GAP_FIELDS,
        usecols=['chrom', 'start', 'end', 'length', 'type', 'bridge'],
        **kwargs)
    if chroms is not None:
        gap = gap[gap.chrom.isin(chroms)]
    return chromsorted(gap)


def _read_bigwig_as_wig(filepath, chrom, start=None, end=None, cachedir=None):
    # https://sebastienvigneau.wordpress.com/2014/01/10/bigwig-to-bedgraph-to-wig/
    # http://redmine.soe.ucsc.edu/forum/index.php?t=msg&goto=5492&S=2925a24be1c20bb064fc09bd054f862d
    cmd = ['bigWigToWig',
           '-chrom={}'.format(chrom)]
    if start is not None:
        cmd += ['-start={}'.format(start)]
    if end is not None:
        cmd += ['-end={}'.format(end)]
    if cachedir is not None:
        cmd += ['-udcDir={}'.format(cachedir)]

    with tempfile.NamedTemporaryFile('w+t') as fh:
        cmd += [filepath, fh.name]
        run(cmd, raises=True)
        fh.flush()

        trackline = fh.readline().split()
        if trackline[0] == '#bedGraph':
            info = {'type': 'bedGraph'}
            out = pd.read_csv(fh, sep='\t', names=['chrom', 'start', 'end', 'value'])
        else:
            tracktype = trackline[0]
            info = dict([kv.split('=') for kv in trackline[1:]])
            info['type'] = tracktype
            for key in ['start', 'step', 'span']:
                if key in info: info[key] = int(info[key])
            if tracktype == 'fixedStep':
                out = pd.read_csv(fh, sep='\t', names=['value'])
            else:
                out = pd.read_csv(fh, sep='\t', names=['start', 'value'])

    return info, out


def read_bigwig_binned(filepath, chrom, start, end, nbins=1, aggfunc=None, cachedir=None):
    """
    Get summary data from bigWig for indicated region, broken into ``nbins`` equal parts.
    
    Parameters
    ----------
    filepath : str
        Path to bigwig file.
    chrom : str
        Chromosome label.
    start, end : int
        Coordinates (zero-based).
    nbins : int, optional
        Number of bins. Default is to summarize the whole region.
    aggfunc : str, optional
        Aggregation method (summary statistic). One of:
        'mean' - average value in region (default)
        'std' - standard deviation in region
        'min' - minimum value in region
        'max' - maximum value in region
        'coverage' - % of region that is covered
    
    """
    cmd = ['bigWigSummary',
           filepath, chrom, str(start+1), str(end), str(nbins)]
    if aggfunc is not None:
        cmd += ['-type={}'.format(aggfunc)]
    if cachedir is not None:
        cmd += ['-udcDir={}'.format(cachedir)]
    out = run(cmd, raises=True)
    return pd.read_csv(StringIO(out), sep='\t', na_values='n/a', header=None).iloc[0].values


def read_bigwig(fp, chrom, start=None, end=None, cachedir=None, as_wiggle=False):
    if as_wiggle:
        return _read_bigwig_as_wig(fp, chrom, start, end, cachedir)

    cmd = [
        'bigWigToBedGraph',
        '-chrom={}'.format(chrom)
    ]
    if start is not None:
        cmd += ['-start={}'.format(start)]
    if end is not None:
        cmd += ['-end={}'.format(end)]
    if cachedir is not None:
        cmd += ['-udcDir={}'.format(cachedir)]

    with tempfile.NamedTemporaryFile('w+t') as fh:
        cmd += [fp, fh.name]
        run(cmd, raises=True)
        fh.flush()
        bg = pd.read_csv(fh, sep='\t', names=['chrom', 'start', 'end', 'value'])

    return bg


def read_tabix(fp, chrom=None, start=None, end=None):
    with closing(pysam.TabixFile(fp)) as f:
        names = list(f.header) or None
        df = pd.read_csv(
            StringIO('\n'.join(f.fetch(chrom, start, end))),
            sep='\t', header=None, names=names)
    return df


def read_bam(fp, chrom=None, start=None, end=None):
    with closing(pysam.AlignmentFile(fp), 'rb') as f:
        bam_iter = f.fetch(chrom, start, end)
        records = [(s.qname, s.flag, s.rname, s.pos, s.mapq,
                    s.cigarstring if s.mapq != 0 else np.nan,
                    s.rnext, s.pnext, s.tlen, s.seq, s.qual,
                    json.dumps(OrderedDict(s.tags))) for s in bam_iter]
        df = pd.DataFrame(records, columns=BAM_FIELDS)
    return df



def read_cytoband(filepath_or_fp):
    return pd.read_csv(
        filepath_or_fp, sep='\t', compression='gzip',
        names=['chrom', 'start', 'end', 'name', 'gieStain'])


def read_fasta(*filepaths, **kwargs):
    records = OrderedDict()
    for filepath in filepaths:
        records.update(pyfaidx.Fasta(filepath, **kwargs))
    return records
