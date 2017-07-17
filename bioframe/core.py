from __future__ import division, print_function
import re

import numpy as np
import pandas as pd


def atoi(s):
    return int(s.replace(',', ''))


def natsort_key(s, _NS_REGEX=re.compile(r'(\d+)', re.U)):
    return tuple([int(x) if x.isdigit() else x for x in _NS_REGEX.split(s) if x])


def natsorted(iterable):
    return sorted(iterable, key=natsort_key)


def argnatsort(array):
    array = np.asarray(array)
    if not len(array): return np.array([], dtype=int)
    cols = tuple(zip(*(natsort_key(x) for x in array)))
    return np.lexsort(cols[::-1])  # numpy's lexsort is ass-backwards


def parse_humanized(s):
    _NUMERIC_RE = re.compile('([0-9,.]+)')
    _, value, unit = _NUMERIC_RE.split(s.replace(',', ''))
    if not len(unit):
        return int(value)
    
    value = float(value)
    unit = unit.upper().strip()
    if unit in ('K', 'KB'):
        value *= 1000
    elif unit in ('M', 'MB'):
        value *= 1000000
    elif unit in ('G', 'GB'):
        value *= 1000000000
    else:
        raise ValueError("Unknown unit '{}'".format(unit))
    return int(value)


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
            ('HYPHEN', r'-'),
            ('COORD',  r'[0-9,]+(\.[0-9]*)?(?:[a-z]+)?'),
            ('OTHER',  r'.+')
        ]
        tok_regex = r'\s*' + r'|\s*'.join(
            r'(?P<%s>%s)' % pair for pair in token_spec)
        tok_regex = re.compile(tok_regex, re.IGNORECASE)
        for match in tok_regex.finditer(s):
            typ = match.lastgroup
            yield typ, match.group(typ)


    def _check_token(typ, token, expected):
        if typ is None:
            raise ValueError('Expected {} token missing'.format(' or '.join(expected)))
        else:
            if typ not in expected:
                raise ValueError('Unexpected token "{}"'.format(token))


    def _expect(tokens):
        typ, token = next(tokens, (None, None))
        _check_token(typ, token, ['COORD'])
        start = parse_humanized(token)

        typ, token = next(tokens, (None, None))
        _check_token(typ, token, ['HYPHEN'])

        typ, token = next(tokens, (None, None))
        if typ is None:
            return start, None

        _check_token(typ, token, ['COORD'])    
        end = parse_humanized(token)
        if end < start:
            raise ValueError('End coordinate less than start')

        return start, end

    parts = s.split(':')
    chrom = parts[0].strip()
    if not len(chrom):
        raise ValueError("Chromosome name cannot be empty")
    if len(parts) < 2:
        return (chrom, None, None)
    start, end = _expect(_tokenize(parts[1]))
    return (chrom, start, end)


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


def bedslice(grouped, chrom, start, end):
    """Assumes no proper nesting of intervals"""
    chromdf = grouped.get_group(chrom)
    lo = chromdf['end'].values.searchsorted(start, side='right')
    hi = lo + chromdf['start'].values[lo:].searchsorted(end, side='left')
    return chromdf.iloc[lo:hi]


def bg2slice_frame(bg2, region1, region2):
    """
    Slice a dataframe with columns ['chrom1', 'start1', 'end1', 'chrom2',
    'start2', 'end2']. Assumes no proper nesting of intervals.
    """
    chrom1, start1, end1 = region1
    chrom2, start2, end2 = region2
    if end1 is None:
        end1 = np.inf
    if end2 is None:
        end2 = np.inf
    out = bg2[(bg2['chrom1'] == chrom1) & 
              (bg2['start1'] >= start1) &
              (bg2['end1'] < end1) &
              (bg2['chrom2'] == chrom2) & 
              (bg2['start2'] >= start2) &
              (bg2['end2'] < end2)]
    return out


def bedslice_frame(bed, region):
    """
    Slice a dataframe with sorted columns ['chrom', 'start', 'end'].
    Assumes no proper nesting of intervals.
    """
    chrom, start, end = region
    grouped = bed.groupby('chrom')
    chromdf = grouped.get_group(chrom)
    lo = chromdf['end'].values.searchsorted(start, side='right')
    if end is not None:
        hi = lo + chromdf['start'].values[lo:].searchsorted(end, side='left')
    else:
        hi = None
    return chromdf.iloc[lo:hi]


def bedslice_series(bed, region):
    """
    Slice a series multi-indexed by ['chrom', 'start', 'end'].
    Assumes no proper nesting of intervals.
    """
    chrom, start, end = region
    return bed.loc[chrom].loc[start:end]


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


# class Genome:
#     '''
#     Tasks: 
    
#     '''
#     def __init__(self,
#             chromsizes='',
#             centromeres=None,
#             bins=None,
#             fasta_path=None,
#             mapped_frac=False,
#             GC=False,
#             ):
#         pass
    
#     def from_cache(path):
#         pass
    
#     @property
#     def chroms(self):
#         return self._chroms
    
#     @property
#     def chromarms(self):
#         return self._chromarms
 
#     @property
#     def bins(self):
#         return self._bins
