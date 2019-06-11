# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
import subprocess
import re
import os

import numpy as np
import pandas as pd

from .region import parse_region


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


def _find_block_span(arr, val):
    '''Find the first and the last occurence + 1 of the value in the array.
    '''
    # it can be done via bisection, but for now BRUTE FORCE
    block_idxs = np.where(arr==val)[0]
    lo, hi = block_idxs[0], block_idxs[-1]+1
    return lo,hi


def bedbisect(bedf, region):
    """Return the span of a block of rows corresponding to
    the genomic region.
    Rows must be sorted by `start` and `end`;
    `chrom` must be grouped, but does not have to be sorted.

    """
    chrom, start, end = parse_region(region)

    lo, hi = _find_block_span(bedf.chrom.values, chrom)

    lo += bedf['end'].values[lo:hi].searchsorted(start, side='right')
    if end is not None:
        hi = lo + bedf['start'].values[lo:hi].searchsorted(end, side='left')
    else:
        hi = None
    return lo, hi


def bedslice(bedf, region):
    """eeturn a block of rows corresponding to the genomic region.
    Rows must be sorted by `start` and `end`;
    `chrom` must be grouped, but does not have to be sorted.
    """
    lo, hi = bedbisect(bedf, region)
    return bedf.iloc[lo:hi]


def bedslice_series(beds, region):
    """
    Slice a series multi-indexed by ['chrom', 'start', 'end'].
    Assumes no proper nesting of intervals.
    """
    chrom, start, end = region
    return beds.loc[chrom].loc[start:end]


def bg2slice(bg2, region1, region2):
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


def expand_regions(df, pad_bp, chromsizes, side='both', inplace=False):
    if not inplace:
        df = df.copy()

    if side == 'both' or side == 'left':
        df.start = np.maximum(0, df.start.values - pad_bp)

    if side == 'both' or side == 'right':
        df.end = np.minimum(df.chrom.apply(chromsizes.__getitem__),
                            df.end+pad_bp)

    return df


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

    tables : sequence of BED-like ``pd.DataFrame``s.
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
        return pd.concat(
            bychrom(lambda c,x:x, df),
            axis=0,
            ignore_index=reset_index)
    else:
        return pd.concat(
            bychrom(lambda c,x:x, df.sort_values(sort_by, **kw)),
            axis=0,
            ignore_index=reset_index)


def split_chromosomes(chromsizes, split_pos, suffixes=('L', 'R')):
    """
    Split chromosomes into chromosome arms

    Parameters
    ----------
    chromsizes : pandas.Series
        Series mapping chromosomes to lengths in bp
    split_pos : pandas.Series
        Series mapping chromosomes to split locations

    Returns
    -------
    4-column BED-like DataFrame (chrom, start, end, name)
    Arm names are chromosome names + L/R suffix.

    """
    index = chromsizes.index.intersection(split_pos.index)
    left_arm = pd.DataFrame(index=index)
    left_arm['chrom'] = left_arm.index
    left_arm['start'] = 0
    left_arm['end'] = split_pos
    left_arm['name'] = left_arm.index + suffixes[0]
    left_arm = left_arm.reset_index(drop=True)

    right_arm = pd.DataFrame(index=index)
    right_arm['chrom'] = right_arm.index
    right_arm['start'] = split_pos
    right_arm['end'] = chromsizes
    right_arm['name'] = right_arm.index + suffixes[1]
    right_arm = right_arm.reset_index(drop=True)

    arms = pd.concat([left_arm, right_arm], axis=0)
    idx = np.lexsort([arms.name, arms.index])
    arms = arms.iloc[idx].reset_index(drop=True)
    return arms


def binnify(chromsizes, binsize, rel_ids=False):
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
        clen = chromsizes[chrom]
        n_bins = int(np.ceil(clen / binsize))
        binedges = np.arange(0, (n_bins+1)) * binsize
        binedges[-1] = clen
        return pd.DataFrame({
                'chrom': [chrom]*n_bins,
                'start': binedges[:-1],
                'end': binedges[1:],
            }, columns=['chrom', 'start', 'end'])

    bintable = pd.concat(map(_each, chromsizes.keys()),
                               axis=0, ignore_index=True)

    if rel_ids:
        bintable['rel_id'] = bintable.groupby('chrom').cumcount()

    # if as_cat:
    #     bintable['chrom'] = pd.Categorical(
    #         bintable['chrom'],
    #         categories=list(chromsizes.keys()),
    #         ordered=True)

    return bintable


def digest(fasta_records, enzyme):
    """
    Divide a genome into restriction fragments.

    Parameters
    ----------
    fasta_records : OrderedDict
        Dictionary of chromosome names to sequence records.
    enzyme: str
        Name of restriction enzyme.

    Returns
    -------
    Dataframe with columns: 'chrom', 'start', 'end'.

    """
    import Bio.Restriction as biorst
    import Bio.Seq as bioseq
    # http://biopython.org/DIST/docs/cookbook/Restriction.html#mozTocId447698
    chroms = fasta_records.keys()
    try:
        cut_finder = getattr(biorst, enzyme).search
    except AttributeError:
        raise ValueError('Unknown enzyme name: {}'.format(enzyme))

    def _each(chrom):
        seq = bioseq.Seq(str(fasta_records[chrom]))
        cuts = np.r_[0, np.array(cut_finder(seq)) + 1, len(seq)].astype(int)
        n_frags = len(cuts) - 1

        frags = pd.DataFrame({
            'chrom': [chrom] * n_frags,
            'start': cuts[:-1],
            'end': cuts[1:]},
            columns=['chrom', 'start', 'end'])
        return frags

    return pd.concat(map(_each, chroms), axis=0, ignore_index=True)


def frac_mapped(bintable, fasta_records):
    def _each(bin):
        s = str(fasta_records[bin.chrom][bin.start:bin.end])
        nbases = len(s)
        n = s.count('N')
        n += s.count('n')
        return (nbases - n) / nbases
    return bintable.apply(_each, axis=1)


def frac_gc(bintable, fasta_records, mapped_only=True):
    def _each(chrom_group):
        chrom = chrom_group.name
        seq = fasta_records[chrom]
        gc = []
        for _, bin in chrom_group.iterrows():
            s = str(seq[bin.start:bin.end])
            g = s.count('G')
            g += s.count('g')
            c = s.count('C')
            c += s.count('c')
            nbases = len(s)
            if mapped_only:
                n = s.count('N')
                n += s.count('n')
                nbases -= n
            gc.append((g + c) / nbases if nbases > 0 else np.nan)
        return gc
    out = bintable.groupby('chrom', sort=False).apply(_each)
    return pd.Series(data=np.concatenate(out), index=bintable.index)


def frac_gene_coverage(bintable, mrna):

    if isinstance(mrna, six.string_types):
        from .resources import fetch_ucsc_mrna
        mrna = fetch_ucsc_mrna(mrna).rename(
            columns={'tName': 'chrom', 'tStart': 'start', 'tEnd': 'end'})

    mrna = mrna.sort_values(['chrom','start','end']).reset_index(drop=True)

    with tsv(bintable) as a, tsv(mrna) as b:
        cov = bedtools.coverage(a=a.name, b=b.name)

    bintable = bintable.copy()
    bintable['gene_count'] = cov.iloc[:,-4]
    bintable['gene_coverage'] = cov.iloc[:,-1]

    return bintable
