from __future__ import division, print_function
import numpy as np
import pandas as pd


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
