from __future__ import division, print_function, unicode_literals
from contextlib import contextmanager, closing
from collections import OrderedDict
import subprocess
import itertools
import functools
import tempfile
import json
import re
import os
import io
from io import StringIO
try:
    from cStringIO import StringIO as BytesIO
except ImportError:
    from io import BytesIO
import six

import numpy as np
import pandas
import pysam

__version__ = '0.1.0'

settings = {
    'bedtools_path': '',
    'kenttools_path': '/net/proteome/home/nezar/local/tools/kenttools/bin',
}


def _is_text_or_bytes(s):
    return isinstance(s, (six.text_type, six.binary_type))


# https://genome.ucsc.edu/FAQ/FAQformat.html
BED12_FIELDS = ['chrom', 'start', 'end',
                'name', 'score', 'strand',
                'thickStart', 'thickEnd', 'rgb',
                'blockCount', 'blockSizes', 'blockStarts']
BED_FIELDS = BED12_FIELDS[:6]

BEDPE_FIELDS = ['chrom1', 'start1', 'end1',
                'chrom2', 'start2', 'end2',
                'name', 'score', 'strand1', 'strand2']
GFF_FIELDS = ['chrom', 'source', 'feature', 'start', 'end',
              'score', 'strand', 'frame', 'attributes']

PGSNP_FIELDS = ['chrom', 'start', 'end', 'name',
                'alleleCount', 'alleleFreq', 'alleleScores']
BEDRNAELEMENTS_FIELDS = ['chrom', 'start', 'end', 'name', 'score', 'strand',
              'level', 'signif', 'score2']
NARROWPEAK_FIELDS = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                     'fc', '-log10p', '-log10q', 'relSummit']
BROADPEAK_FIELDS = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                     'fc', '-log10p', '-log10q']
GAPPEDPEAK_FIELDS = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                     'thickStart', 'thickEnd', 'rgb',
                     'blockCount', 'blockSizes', 'blockStarts',
                     'fc', '-log10p', '-log10q']
GAP_FIELDS = ['bin', 
              'chrom', 'start', 'end', 
              'ix', 'n', 
              'length', 'type', 'bridge']

# http://ga4gh.org/#/fileformats-team
BAM_FIELDS = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR',
              'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'TAGs']
VCF_FIELDS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

SCHEMAS = {
    'bed': BED_FIELDS,
    'bedGraph': BED_FIELDS[:3],
    'bed3': BED_FIELDS[:3],
    'bed4': BED_FIELDS[:4],
    'bed5': BED_FIELDS[:5],
    'bed6': BED_FIELDS,
    'bed9': BED12_FIELDS[:9],
    'bed12': BED12_FIELDS,
    'gff': GFF_FIELDS,
    'gtf': GFF_FIELDS,
    'bedRnaElements': BEDRNAELEMENTS_FIELDS,
    'narrowPeak': NARROWPEAK_FIELDS,
    'broadPeak': BROADPEAK_FIELDS,
    'gappedPeak': GAPPEDPEAK_FIELDS,
    'bam': BAM_FIELDS,
    'vcf': VCF_FIELDS,
}

CHROMINFO_URLS = {
    'hg38': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz',
    'hg19': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz',
    'mm10': 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromInfo.txt.gz',
    'mm9': 'http://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/chromInfo.txt.gz',
}



def atoi(s):
    return int(s.replace(',', ''))


def parse_region_string(s):
    """
    Parse a UCSC-style genomic region string into a triple.
    Parameters
    ----------
    s : str
        UCSC-style string, e.g. "chr5:10,100,000-30,000,000". Ensembl and FASTA
        style sequence names are allowed. End coordinate must be greater than or
        equal to start.
    
    Returns
    -------
    (str, int or None, int or None)
    """
    def _tokenize(s):
        token_spec = [
            ('INT',    r'[0-9,]+'),
            ('ALNUM',  r'[a-zA-z0-9_|]+'),
            ('COLON',  r':'),
            ('HYPHEN', r'-'),
        ]
        tok_regex = r'\s*' + r'|\s*'.join(
            r'(?P<%s>%s)' % pair for pair in token_spec)
        tok_regex = re.compile(tok_regex)
        for match in tok_regex.finditer(s):
            typ = match.lastgroup
            yield typ, match.group(typ)

    def _check_next(tokens, expected):
        try:
            token = next(tokens)
        except StopIteration:
            raise ValueError('Expected {} token missing'.format(expected))
        else:
            if token[0] not in expected:
                raise ValueError('Unexpected token "{}"'.format(token[1]))
        return token[1]

    def _expect(tokens):
        chrom = _check_next(tokens, ['ALNUM', 'INT'])
        try:
            token = next(tokens)
        except StopIteration:
            return (chrom, None, None)
        if token[0] != 'COLON':
            raise ValueError('Got "{}" after chromosome label'.format(token[1]))

        start = atoi(_check_next(tokens, ['INT']))
        _check_next(tokens, ['HYPHEN'])
        end = atoi(_check_next(tokens, ['INT']))
        if end < start:
            raise ValueError('End coordinate less than start')

        return chrom, start, end

    return _expect(_tokenize(s))


def parse_region(reg, chromsizes=None):
    """
    Genomic regions are represented as half-open intervals (0-based starts,
    1-based ends) along the length coordinate of a contig/scaffold/chromosome.
    Parameters
    ----------
    reg : str or tuple
        UCSC-style genomic region string, or
        Triple (chrom, start, end), where ``start`` or ``end`` may be ``None``.
    chromsizes : mapping, optional
        Lookup table of scaffold lengths to check against ``chrom`` and the
        ``end`` coordinate. Required if ``end`` is not supplied.
    Returns
    -------
    A well-formed genomic region triple (str, int, int)
    """
    if isinstance(reg, six.string_types):
        chrom, start, end = parse_region_string(reg)
    else:
        chrom, start, end = reg
        start = int(start) if start is not None else start
        end = int(end) if end is not None else end

    try:
        clen = chromsizes[chrom] if chromsizes is not None else None
    except KeyError:
        raise ValueError("Unknown sequence label: {}".format(chrom))
    
    start = 0 if start is None else start
    if end is None:
        if clen is None:  # TODO --- remove?
            raise ValueError("Cannot determine end coordinate.")
        end = clen

    if end < start:
        raise ValueError("End cannot be less than start")
    
    if start < 0 or (clen is not None and end > clen):
        raise ValueError(
            "Genomic region out of bounds: [{}, {})".format(start, end))
    
    return chrom, start, end



def natsort_key(s, _NS_REGEX=re.compile(r'(\d+)', re.U)):
    return tuple([int(x) if x.isdigit() else x for x in _NS_REGEX.split(s) if x])


def natsorted(iterable):
    return sorted(iterable, key=natsort_key)


def argnatsort(array):
    array = np.asarray(array)
    if not len(array): return np.array([], dtype=int)
    cols = tuple(zip(*(natsort_key(x) for x in array)))
    return np.lexsort(cols[::-1])  # numpy's lexsort is ass-backwards


@contextmanager
def open_stream(fp, mode='rb', gzipped=False, *args, **kwargs):
    """
    Context manager like ``open`` but allows already open file handles and
    file-like objects to pass through without being closed.

    Parameters
    ----------
    fp : string or file-like

    Other parameters
    ----------------
    args, kwargs : tuple, dict
        Extra arguments are passed on to the ``open`` builtin when fp is a
        string.

    """
    own_fh = False
    if _is_text_or_bytes(fp):
        own_fh = True
        if fp.endswith('.gz'): gzipped = True
        fh = open(fp, mode, *args, **kwargs)
    else:
        fh = fp

    try:
        if gzipped:
            if 'r' in mode:
                fz = io.BufferedReader(gzip.GzipFile(fh.name, mode, fileobj=fh))
            else:
                fz = io.BufferedWriter(gzip.GzipFile(fh.name, mode, fileobj=fh))
            yield fz
        else:
            yield fh
    finally:
        if gzipped: # Calling a GzipFile object's close() method does not close fileobj
            fz.close()
        if own_fh:
            fh.close()


def run(cmd, input=None, raises=True, print_cmd=False, max_msg_len=1000):
    if print_cmd:
        print(subprocess.list2cmdline(cmd))

    if input is not None:
        p = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = p.communicate(input.encode('utf-8'))
    else:
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        out, err = p.communicate()


    if raises and p.returncode != 0:
        if len(out) > max_msg_len:
            out = out[:max_msg_len] + b'... [truncated]'
        raise OSError("process failed: %d\n%s\n%s" % (p.returncode,  out.decode('utf-8'), err.decode('utf-8')))

    return out.decode('utf-8')


def read_chromsizes(filepath_or,
                   name_patterns=(r'^chr[0-9]+$', r'^chr[XY]$', r'^chrM$'),
                   all_names=False,
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
    * NCBI assembly terminology: <https://www.ncbi.nlm.nih.gov/grc/help/definitions

    """
    if isinstance(filepath_or, six.string_types) and filepath_or.endswith('.gz'):
        kwargs.setdefault('compression', 'gzip')
    chromtable = pd.read_csv(
        filepath_or, sep='\t', usecols=[0, 1],
        names=['name', 'length'], dtype={'name':str}, **kwargs)
    if not all_names:
        parts = []
        for pattern in name_patterns:
            part = chromtable[chromtable['name'].str.contains(pattern)]
            part = part.iloc[argnatsort(part['name'])]
            parts.append(part)
        chromtable = pd.concat(parts, axis=0)
    chromtable.index = chromtable['name'].values
    return chromtable['length']


def fetch_chromsizes(db, **kwargs):
    """
    Download chromosome sizes from UCSC as a ``pandas.Series``, indexed by
    chromosome label.

    """
    return read_chromsizes(
        'http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/chromInfo.txt.gz'.format(db), 
        **kwargs)


def read_gapfile(filepath_or_fp, chroms=None, **kwargs):
    gap = pandas.read_csv(
        filepath_or_fp,
        sep='\t',
        names=GAP_FIELDS,
        usecols=['chrom', 'start', 'end', 'length', 'type', 'bridge'],
        **kwargs)
    if chroms is not None:
        gap = gap[gap.chrom.isin(chroms)]
    return chromsorted(gap)


def fetch_gaps(db, **kwargs):
    return read_gapfile('http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/gap.txt.gz'.format(db), **kwargs)


def read_gapfile(filepath_or_fp, chroms=None, **kwargs):
    gap = pandas.read_csv(
        filepath_or_fp,
        sep='\t',
        names=GAP_FIELDS,
        usecols=['chrom', 'start', 'end', 'length', 'type', 'bridge'],
        **kwargs)
    if chroms is not None:
        gap = gap[gap.chrom.isin(chroms)]
    return chromsorted(gap)


def read_table(filepath_or, schema=None, **kwargs):
    kwargs.setdefault('sep', '\t')
    kwargs.setdefault('header', None)
    if _is_text_or_bytes(filepath_or) and filepath_or.endswith('.gz'):
        kwargs.setdefault('compression', 'gzip')
    if schema is not None:
        try:
            kwargs.setdefault('names', SCHEMAS[schema])
        except (KeyError, TypeError):
            if isinstance(schema, six.string_types):
                raise ValueError("TSV schema not found: '{}'".format(schema))
            kwargs.setdefault('names', schema)
    return pandas.read_csv(filepath_or, **kwargs)


def _read_bigwig_as_wig(filepath, chrom, start=None, end=None, cachedir=None):
    # https://sebastienvigneau.wordpress.com/2014/01/10/bigwig-to-bedgraph-to-wig/
    # http://redmine.soe.ucsc.edu/forum/index.php?t=msg&goto=5492&S=2925a24be1c20bb064fc09bd054f862d
    cmd = [os.path.join(settings['kenttools_path'], 'bigWigToWig'),
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
            out = pandas.read_csv(fh, sep='\t', names=['chrom', 'start', 'end', 'value'])
        else:
            tracktype = trackline[0]
            info = dict([kv.split('=') for kv in trackline[1:]])
            info['type'] = tracktype
            for key in ['start', 'step', 'span']:
                if key in info: info[key] = int(info[key])
            if tracktype == 'fixedStep':
                out = pandas.read_csv(fh, sep='\t', names=['value'])
            else:
                out = pandas.read_csv(fh, sep='\t', names=['start', 'value'])

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
    cmd = [os.path.join(settings['kenttools_path'], 'bigWigSummary'),
           filepath, chrom, str(start+1), str(end), str(nbins)]
    if aggfunc is not None:
        cmd += ['-type={}'.format(aggfunc)]
    if cachedir is not None:
        cmd += ['-udcDir={}'.format(cachedir)]
    out = run(cmd, raises=True)
    return pandas.read_csv(StringIO(out), sep='\t', na_values='n/a', header=None).iloc[0].values


def read_bigwig(fp, chrom, start=None, end=None, cachedir=None, as_wiggle=False):
    if as_wiggle:
        return _read_bigwig_as_wig(fp, chrom, start, end, cachedir)

    cmd = [
        os.path.join(settings['kenttools_path'], 'bigWigToBedGraph'),
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
        bg = pandas.read_csv(fh, sep='\t', names=['chrom', 'start', 'end', 'value'])

    return bg


def read_tabix(fp, chrom=None, start=None, end=None):
    with closing(pysam.TabixFile(fp)) as f:
        names = list(f.header) or None
        df = pandas.read_csv(
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
        df = pandas.DataFrame(records, columns=BAM_FIELDS)
    return df


def binnify(chromsizes, binsize):
    """
    Divide a genome into evenly sized bins.

    Parameters
    ----------
    chromsizes : Series
        pandas Series indexed by chromosome name with chromosome lengths in bp.
    binsize : int
        size of bins in bp

    Returns
    -------
    Data frame with columns: 'chrom', 'start', 'end'.

    """
    def _each(chrom):
        clen = chromsizes.at[chrom]
        n_bins = int(np.ceil(clen / binsize))
        binedges = np.arange(0, (n_bins+1)) * binsize
        binedges[-1] = clen
        return pandas.DataFrame({
                'chrom': [chrom]*n_bins,
                'start': binedges[:-1],
                'end': binedges[1:],
            }, columns=['chrom', 'start', 'end'])
    return pandas.concat(map(_each, chromsizes.index),
                         axis=0, ignore_index=True)


def bychrom(func, *tables, **kwargs):
    """
    Split one or more bed-like dataframes by chromosome.
    Apply ``func(chrom, *slices)`` to each chromosome slice.
    Yield results.

    Parameters
    ----------
    func : function to apply to split dataframes.
        The expected signature is ``func(chrom, df1[, df2[, ...])``,
        where ``df1, df2, ...`` are subsets of the input dataframes.
        The function can return anything.

    tables : sequence of BED-like ``pandas.DataFrame``s.
        The first column of each dataframe must be chromosome labels,
        unless specified by ``chrom_field``.

    chroms : sequence of str, optional
        Select which chromosome subsets of the data to apply the function to.
        Defaults to all unique chromosome labels in the first dataframe input,
        in natural sorted order.

    chrom_field: str, optional
        Name of column containing chromosome labels.

    ret_chrom : bool, optional (default: False)
        Yield "chromosome, value" pairs as output instead of only values.

    map : callable, optional (default: ``itertools.imap`` or ``map`` in Python 3)
        Map implementation to use.

    Returns
    -------
    Iterator or future that yields the output of running `func` on
    each chromosome

    """
    chroms = kwargs.pop('chroms', None)
    parallel = kwargs.pop('parallel', False)
    ret_chrom = kwargs.pop('ret_chrom', False)
    map_impl = kwargs.pop('map', six.moves.map)

    first = tables[0]
    chrom_field = kwargs.pop('chrom_field', first.columns[0])
    if chroms is None:
        chroms = natsorted(first[chrom_field].unique())

    grouped_tables = [table.groupby(chrom_field) for table in tables]

    def iter_partials():
        for chrom in chroms:
            partials = []
            for gby in grouped_tables:
                try:
                    partials.append(gby.get_group(chrom))
                except KeyError:
                    partials.append(gby.head()[0:0])
            yield partials

    if ret_chrom:
        def run_job(chrom, partials):
            return chrom, func(chrom, *partials)
    else:
        def run_job(chrom, partials):
            return func(chrom, *partials)

    return map_impl(run_job, chroms, iter_partials())


def chromsorted(df, sort_by=None, reset_index=True, **kw):
    """
    Sort bed-like dataframe by chromosome label in "natural" alphanumeric order,
    followed by any columns specified in ``sort_by``.

    """
    if sort_by is None:
        return pandas.concat(
            bychrom(lambda c,x:x, df),
            axis=0,
            ignore_index=reset_index)
    else:
        return pandas.concat(
            bychrom(lambda c,x:x, df.sort_values(sort_by, **kw)),
            axis=0,
            ignore_index=reset_index)


def tsv(df, **kwargs):
    """
    Write ``pandas.DataFrame`` to a temporary tab-delimited file.
    Works in a ``with`` block (file is deleted at context teardown).

    >>> with tsv(df1) as f1, tsv(df2) as f2:
    ...    # something that requires tsv file input (use f or f.name)

    """
    fh = tempfile.NamedTemporaryFile(mode='w+t')
    df.to_csv(fh, sep=str('\t'), index=False, header=False, **kwargs)
    fh.flush()  # DON'T FORGET TO FLUSH!!!
    fh.seek(0)
    return fh


def to_dataframe(text, columns=None):
    # To convert decoded stdout into a dataframe
    return pandas.read_csv(StringIO(text), sep='\t', header=None, names=columns)


def cmd_exists(cmd):
    return any(os.access(os.path.join(path, cmd), os.X_OK)
               for path in os.environ["PATH"].split(os.pathsep))


if cmd_exists("bedtools"):

    def _register(name):
        # Wrapper calls bedtools
        def wrapper(**kwargs):
            columns = kwargs.pop('_schema', None)
            run_kws = {kw[1:]: kwargs.pop(kw) for kw in list(kwargs.keys())
                           if kw.startswith('_')}

            cmd = [os.path.join(settings['bedtools_path'], 'bedtools'), name]
            for k, v in kwargs.items():
                if isinstance(v, bool):
                    if not v: continue
                    cmd.append('-{}'.format(k))
                else:
                    cmd.append('-{}'.format(k))
                    cmd.append(str(v))
            out = run(cmd, **run_kws)
            return to_dataframe(out, columns=columns)

        # Call once to generate docstring from usage text
        p = subprocess.Popen(
                [os.path.join(settings['bedtools_path'],'bedtools'), name, '-h'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        _, err = p.communicate()
        wrapper.__doc__ = err.decode('utf-8')
        wrapper.__name__ = str(name)

        return staticmethod(wrapper)


    class bedtools(object):
        intersect = _register('intersect')
        window = _register('window')
        closest = _register('closest')
        coverage = _register('coverage')
        map = _register('map')
        genomecov = _register('genomecov')
        merge = _register('merge')
        cluster = _register('cluster')
        complement = _register('complement')
        subtract = _register('subtract')
        slop = _register('slop')
        flank = _register('flank')
        sort = _register('sort')
        random = _register('random')
        shuffle = _register('shuffle')
        annotate = _register('annotate')
        jaccard = _register('jaccard')


class IndexedBedLike(object):
    """BED intersection using pandas"""
    def __init__(self, bed):
        # create two sorted lookup indexes
        self.lookup_head = bed.set_index(
            ['chrom',   'end'], verify_integrity=True).sortlevel()
        self.lookup_tail = bed.set_index(
            ['chrom', 'start'], verify_integrity=True).sortlevel()

    def intersect(self, qchrom, qstart, qend):
        # fetch all intervals that terminate inside the query window, no matter where they begin
        head = self.lookup_head[(qchrom, qstart):(qchrom, qend)].reset_index()

        # fetch all intervals that begin inside the query window...
        tail = self.lookup_tail[(qchrom, qstart):(qchrom, qend)].reset_index()
        # ...and terminate outside it
        tail = tail[tail['end'] > qend]

        return pandas.concat((head, tail), axis=0)


def bedslice(grouped, chrom, start, end):
    """Assumes no proper nesting of intervals"""
    chromdf = grouped.get_group(chrom)
    lo = chromdf['end'].values.searchsorted(start, side='right')
    hi = lo + chromdf['start'].values[lo:].searchsorted(end, side='left')
    return chromdf.iloc[lo:hi]


def frac_mapped(bintable, fasta_records):
    def _each(bin):
        s = str(fasta_records[bin.chrom][bin.start:bin.end])
        nbases = len(s)
        n = s.count('N')
        n += s.count('n')
        return (nbases - n) / nbases
    return bintable.apply(_each, axis=1)


def frac_gc(bintable, fasta_records, mapped_only=True):
    def _each(bin):
        s = str(fasta_records[bin.chrom][bin.start:bin.end])
        g = s.count('G')
        g += s.count('g')
        c = s.count('C')
        c += s.count('c')
        nbases = len(s)
        if mapped_only:
            n = s.count('N')
            n += s.count('n')
            nbases -= n
        return (g + c) / nbases if nbases > 0 else np.nan
    return bintable.apply(_each, axis=1)
