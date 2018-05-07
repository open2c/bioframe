from __future__ import division, print_function
from collections import OrderedDict, Mapping
import subprocess
import tempfile
import json
import six
import os
import io

import numpy as np
import pandas as pd
import pyfaidx
import pysam
from .process import run
from ..core import argnatsort
from ..schemas import SCHEMAS, GAP_FIELDS, UCSC_MRNA_FIELDS


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
    return gap


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
    return pd.read_csv(io.StringIO(out), sep='\t', na_values='n/a', header=None).iloc[0].values


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
            io.StringIO('\n'.join(f.fetch(chrom, start, end))),
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


def read_ucsc_mrnafile(filepath_or_fp, chroms=None, **kwargs):
    mrna = pd.read_csv(
        filepath_or_fp,
        sep='\t',
        names=UCSC_MRNA_FIELDS,
        #usecols=['chrom', 'start', 'end', 'length', 'type', 'bridge'],
        **kwargs)
    if chroms is not None:
        mrna = mrna[mrna.chrom.isin(chroms)]
    return mrna


def load_fasta(filepath_or, engine='pysam', **kwargs):
    """
    Load lazy fasta sequences from an indexed fasta file (optionally compressed)
    or from a collection of uncompressed fasta files.
    
    Parameters
    ----------
    filepath_or : str or iterable
        If a string, a filepath to a single `.fa` or `.fa.gz` file. Assumed to 
        be accompanied by a `.fai` index file. Depending on the engine, the 
        index may be created on the fly, and some compression formats may not 
        be supported. If not a string, an iterable of fasta file paths each 
        assumed to contain a single sequence.
    engine : {'pysam', 'pyfaidx'}, optional
        Module to use for loading sequences.
    kwargs : optional
        Options to pass to ``pysam.FastaFile`` or ``pyfaidx.Fasta``.
    
    Returns
    -------
    OrderedDict of (lazy) fasta records.
    
    Notes
    -----
    * pysam/samtools can read .fai and .gzi indexed files, I think.
    * pyfaidx can handle uncompressed and bgzf compressed files.
    
    """
    class PysamFastaRecord(object):
        def __init__(self, ff, ref):
            self.ff = ff
            if ref not in ff.references:
                raise KeyError(
                    "Reference name '{}' not found in '{}'".format(ref, ff))
            self.ref = ref
        def __getitem__(self, key):
            if isinstance(key, slice):
                start, stop = key.start, key.stop
            else:
                start = key
                stop = key + 1
            return self.ff.fetch(self.ref, start, stop)

    is_multifile = not isinstance(filepath_or, six.string_types)
    records = OrderedDict()
    
    if engine == 'pysam':
        try:
            import pysam
        except ImportError:
            raise ImportError("pysam is required to use engine='pysam'")
        
        if is_multifile:
            for onefile in filepath_or:
                ff = pysam.FastaFile(onefile, **kwargs)
                name = ff.references[0]
                records[name] = PysamFastaRecord(ff, name) 
        else:
            ff = pysam.FastaFile(filepath_or, **kwargs)
            for name in ff.references:
                records[name] = PysamFastaRecord(ff, name) 

    elif engine == 'pyfaidx':
        try:
            import pyfaidx
        except ImportError:
            raise ImportError("pyfaidx is required to use engine='pyfaidx'")

        if is_multifile:
            for onefile in filepath_or:
                ff = pyfaidx.Fasta(onefile, **kwargs)
                name = next(iter(ff.keys()))
                records[name] = ff[name]
        else:
            ff = pyfaidx.Fasta(filepath_or, **kwargs)
            for name in ff.keys():
                records[name] = ff[name]
    
    else:
        raise ValueError("engine must be 'pysam' or 'pyfaidx'")
    
    return records


def to_bigwig(df, chromsizes, outpath, value_field=None):
    """
    Save a bedGraph-like dataframe as a binary BigWig track.

    Parameters
    ----------
    df : pandas.DataFrame
        Data frame with columns 'chrom', 'start', 'end' and one or more value 
        columns
    chromsizes : pandas.Series
        Series indexed by chromosome name mapping to their lengths in bp
    outpath : str
        The output BigWig file path
    value_field : str, optional
        Select the column label of the data frame to generate the track. Default
        is to use the fourth column.

    """
    is_bedgraph = True
    for col in ['chrom', 'start', 'end']:
        if col not in df.columns:
            is_bedgraph = False
    if len(df.columns) < 4:
        is_bedgraph = False

    if not is_bedgraph:
        raise ValueError(
            "A bedGraph-like DataFrame is required, got {}".format(
                df.columns))

    if value_field is None:
        value_field = df.columns[3]

    columns = ['chrom', 'start', 'end', value_field]
    bg = df[columns].copy()
    bg['chrom'] = bg['chrom'].astype(str)
    bg = bg.sort_values(['chrom', 'start', 'end'])

    with tempfile.NamedTemporaryFile(suffix='.bg') as f, \
         tempfile.NamedTemporaryFile('wt', suffix='.chrom.sizes') as cs:

        chromsizes.to_csv(cs, sep='\t')
        cs.flush()

        bg.to_csv(
            f.name,
            sep='\t', 
            columns=columns,
            index=False,
            header=False,
            na_rep='nan')

        run(['bedGraphToBigWig', f.name, cs.name, outpath], 
            print_cmd=True)


def to_bigbed(df, chromsizes, outpath, schema='bed6'):
    """
    Save a bedGraph-like dataframe as a binary BigWig track.

    Parameters
    ----------
    df : pandas.DataFrame
        Data frame with columns 'chrom', 'start', 'end' and one or more value 
        columns
    chromsizes : pandas.Series
        Series indexed by chromosome name mapping to their lengths in bp
    outpath : str
        The output BigWig file path
    value_field : str, optional
        Select the column label of the data frame to generate the track. Default
        is to use the fourth column.

    """
    import tempfile
    import subprocess
    is_bed6 = True
    for col in ['chrom', 'start', 'end', 'name', 'score', 'strand']:
        if col not in df.columns:
            is_bed6 = False
    if len(df.columns) < 6:
        is_bed6 = False

    if not is_bed6:
        raise ValueError(
            "A bed6-like DataFrame is required, got {}".format(
                df.columns))

    columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    bed = df[columns].copy()
    bed['chrom'] = bed['chrom'].astype(str)
    bed = bed.sort_values(['chrom', 'start', 'end'])

    with tempfile.NamedTemporaryFile(suffix='.bed') as f, \
         tempfile.NamedTemporaryFile('wt', suffix='.chrom.sizes') as cs:

        chromsizes.to_csv(cs, sep='\t')
        cs.flush()

        bed.to_csv(
            f.name,
            sep='\t', 
            columns=columns,
            index=False,
            header=False,
            na_rep='nan')

        p = subprocess.run([
                'bedToBigBed', 
                '-type={}'.format(schema), 
                f.name, 
                cs.name, 
                outpath
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    return p


def to_parquet(pieces, outpath, row_group_size=None, compression='snappy',
               use_dictionary=True, version=2.0, **kwargs):
    """
    Save an iterable of dataframe chunks to a single Apache Parquet file. For 
    more info about Parquet, see https://arrow.apache.org/docs/python/parquet.html.

    Parameters
    ----------
    pieces : DataFrame or iterable of DataFrame
        Chunks to write
    outpath : str
        Path to output file
    row_group_size : int
        Number of rows per row group
    compression : {'snappy', 'gzip', 'brotli', 'none'}, optional
        Compression algorithm. Can be set on a per-column basis with a 
        dictionary of column names to compression lib.
    use_dictionary : bool, optional
        Use dictionary encoding. Can be set on a per-column basis with a list
        of column names.

    See also
    --------
    pyarrow.parquet.write_table
    pyarrow.parquet.ParquetFile
    fastparquet

    """
    try:
        import pyarrow.parquet
        import pyarrow as pa
    except ImportError:
        raise ImportError('Saving to parquet requires the `pyarrow` package')

    if isinstance(pieces, pd.DataFrame):
        pieces = (pieces,)

    try:
        for i, piece in enumerate(pieces):
            table = pa.Table.from_pandas(piece, preserve_index=False)
            if i == 0:
                writer = pa.parquet.ParquetWriter(
                    outpath,
                    table.schema,
                    compression=compression,
                    use_dictionary=use_dictionary,
                    version=version,
                    **kwargs)
            writer.write_table(table, row_group_size=row_group_size)
    finally:
        writer.close()
