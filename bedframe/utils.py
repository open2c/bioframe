from __future__ import division, print_function, unicode_literals

import numpy as np
import pandas
import re

import six

def _is_text_or_bytes(s):
    return isinstance(s, (six.text_type, six.binary_type))

def natsort_key(s, _NS_REGEX=re.compile(r'(\d+)', re.U)):
    return tuple([int(x) if x.isdigit() else x for x in _NS_REGEX.split(s) if x])

def natsorted(iterable):
    return sorted(iterable, key=natsort_key)

def argnatsort(array):
    array = np.asarray(array)
    if not len(array): return np.array([], dtype=int)
    cols = tuple(zip(*(natsort_key(x) for x in array)))
    return np.lexsort(cols[::-1])  # numpy's lexsort is ass-backwards

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

