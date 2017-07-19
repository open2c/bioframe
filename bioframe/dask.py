from __future__ import division, print_function, absolute_import
from collections import OrderedDict

import numpy as np
import pandas as pd
import numba

import pypairix
import pysam
from dask.base import tokenize
import dask.dataframe as dd
import dask.array as da
import dask


def bin2start(k):
    lev = np.floor(np.log2(7*k + 1)/3).astype(int)
    sl = 2**(29 - 3*lev)
    ol = (2**(3*lev) - 1)//7
    start = (k - ol) * sl
    end = (k - ol+1) * sl
    return start

LEVEL = {}
LEVEL[0] = bin2start(np.arange(1, 9))
LEVEL[1] = bin2start(np.arange(9, 73))
LEVEL[2] = bin2start(np.arange(73,585))
LEVEL[3] = bin2start(np.arange(585,4681))
LEVEL[4] = bin2start(np.arange(4681,37449))


@numba.jit("int32(int32, int32)")
def reg2bin(beg, end):
    end -= 1
    if beg >> 14 == end >> 14: 
        return ((1 << 15)-1) // 7 + (beg >> 14)
    if beg >> 17 == end >> 17: 
        return ((1 << 12)-1) // 7 + (beg >> 17)
    if beg >> 20 == end >> 20: 
        return ((1 << 9)-1) // 7 + (beg >> 20)
    if beg >> 23 == end >> 23: 
        return ((1 << 6)-1) // 7 + (beg >> 23)
    if beg >> 26 == end >> 26: 
        return ((1 << 3)-1) // 7 + (beg >> 26)
    return 0


@numba.jit("int32(int32, int32)")
def reg2bins(rbeg, rend):
    
    MAX_BIN = ((1 << 18) - 1) // 7
    lst = []

    rend -= 1    
    k = 1 + (rbeg >> 26)
    while k <= (1 + (rend >> 26)):
        k += 1
        lst.append(k)
    
    k = 9 + (rbeg >> 23)
    while k <= (9 + (rend >> 23)):
        k += 1
        lst.append(k)
        
    k = 73 + (rbeg >> 20)
    while k <= (73 + (rend >> 20)):
        k += 1
        lst.append(k)
    
    k = 585 + (rbeg >> 17)
    while k <= (585 + (rend >> 17)):
        k += 1
        lst.append(k)
    
    k = 4681 + (rbeg >> 14)
    while k <= (4681 + (rend >> 14)):
        k += 1
        lst.append(k)
    
    return lst


def range_partition(start, stop, step):
    return ((i, min(i+step, stop))
                for i in range(start, stop, step))


def _fetch_region(filepath, chromsizes, slc, chrom1, chrom2=None, 
                 columns=None, usecols=None, meta=None):
    if chrom2 is None:
        chrom2 = chrom1
    if slc is None:
        start, end = 0, chromsizes[chrom1]
    else:
        start, end = slc.start, slc.stop
    f = pypairix.open(filepath, 'r')
    df = pd.DataFrame.from_records(
        f.query2D(chrom1, start, end, chrom2, 0, chromsizes[chrom2]), 
        columns=columns)

    if not len(df):
        df = meta.copy()
    elif usecols is not None:
        usecols = set(usecols)
        df = df[[col for col in meta.columns if col in usecols]]
    
    for col, dt in meta.dtypes.items():
        df.loc[:, col] = df.loc[:, col].astype(dt)
    
    # nasty hack!
    if len(df) == 0:
        class fake_loc:
            def __init__(self, obj):
                self.obj = obj
            def __call__(self, *args):
                return self.obj
            def __getitem__(self, *args):
                return self.obj
        df._loc = fake_loc(df)
    return df


def daskify_pairix_block(filepath, chromsizes, chrom1, chrom2=None, 
                         columns=None, dtypes=None, usecols=None, 
                         chunk_level=2):
    nrows = chromsizes[chrom1]
    meta = pd.read_csv(
        filepath, 
        sep='\t', 
        comment='#', 
        header=None,
        names=columns,
        dtype=dtypes,
        usecols=usecols,
        iterator=True).read(1024).iloc[0:0]

    # Make a unique task name
    token = tokenize(filepath, chromsizes, chrom1, chrom2, 
                     columns, dtypes, usecols, chunk_level)
    task_name = 'daskify-pairix-' + token

    # Build the task graph
    divisions = []
    dsk = {}
    edges = LEVEL[chunk_level]
    edges = edges[:np.searchsorted(edges, nrows)]
    if edges[-1] != nrows:
        edges = np.r_[edges, nrows]
    spans = zip(edges[:-1], edges[1:])
    for i, (lo, hi) in enumerate(spans):
        divisions.append(hi-1)
        slc = slice(lo, hi)
        dsk[task_name, i] = (_fetch_region, 
                             filepath, chromsizes, slc, 
                             chrom1, chrom2, columns, usecols, meta)
        
    # Generate ddf from dask graph
    return dd.DataFrame(dsk, task_name, meta, tuple(divisions))


def daskify_pairix(filepath, chromsizes, **kwargs):
    f = pypairix.open(filepath)
    blocks = [s.split('|') for s in f.get_blocknames()]
    d = OrderedDict()
    for chrom1, chrom2 in blocks:
        if chrom1 in chromsizes and chrom2 in chromsizes:
            d[chrom1, chrom2] = daskify_pairix_block(
                filepath, chromsizes, chrom1, chrom2, **kwargs)
    return d
