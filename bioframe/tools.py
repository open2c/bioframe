# -*- coding: utf-8 -*-
from __future__ import division, print_function
import subprocess
import os

import numpy as np
import pandas as pd

from .io.process import cmd_exists, run, to_dataframe, tsv
from .io.resources import fetch_ucsc_mrna



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
            s = seq[bin.start:bin.end]
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


def frac_gene_coverage(bintable, db):
    mrna = fetch_ucsc_mrna(db)
    mrna = mrna[['tName','tStart','tEnd']]
    mrna.columns=['chrom','start','end']
    mrna = mrna.sort_values(['chrom','start','end']).reset_index(drop=True)

    with tsv(bintable) as a, tsv(mrna) as b:
        cov = bedtools.coverage(a=a.name, b=b.name)

    bintable = bintable.copy()
    bintable['gene_count'] = cov.iloc[:,-4]
    bintable['gene_coverage'] = cov.iloc[:,-1]
                        
    return bintable



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


def parse_gtf_attributes(attrs, kv_sep='=', item_sep=';', **kwargs):
    item_lists = attrs.str.split(item_sep)
    item_lists = item_lists.apply(
        lambda items: [item.strip().split(kv_sep) for item in items]
    )
    item_lists = item_lists.apply(
        lambda items: [map(str.strip, item) for item in items if len(item) == 2]
    )
    kv_records = item_lists.apply(dict)
    return pd.DataFrame.from_records(kv_records, **kwargs)

# def bedslice(df, chrom, start, end, is_sorted=True, has_overlaps=False, 
#     allow_partial=False, trim_partial=False):
#     '''
#     Extract all bins of `chrom` between `start` and `end`.
    
# #     '''
# #     pass

# # def downsample():
# #     '''
# #     '''

# # def upsample():
# #     '''
# #     '''

# # def rle():
# #     '''
# #     run-length encode a bed chunk
# #     '''

# # def pile():
# #     '''
# #     make a pile up
#     '''

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

        return pd.concat((head, tail), axis=0)



if cmd_exists("bedtools"):

    def _register(name):
        # Wrapper calls bedtools
        def wrapper(**kwargs):
            columns = kwargs.pop('_schema', None)
            run_kws = {kw[1:]: kwargs.pop(kw) for kw in list(kwargs.keys())
                           if kw.startswith('_')}

            cmd = ['bedtools', name]
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
                ['bedtools', name, '-h'],
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


def intersect(bed1, bed2, overlap=True, outer_join=False, v=False, sort=False, suffixes=('_x', '_y')):
    
    # hacky, but we don't want to use suffixes when using -v mode
    if v:
        suffixes = ('',)
    
    bed1_extra = bed1[bed1.columns.difference(['chrom', 'start', 'end'])]
    bed2_extra = bed2[bed2.columns.difference(['chrom', 'start', 'end'])]

    left = bed1[['chrom', 'start', 'end']].copy()
    left['index'] = left.index

    right = bed2[['chrom', 'start', 'end']].copy()
    right['index'] = right.index

    bt_kwargs = {
        'v': v,
        'nonamecheck': False,
    }

    if outer_join:
        if overlap:
            bt_kwargs['wao'] = True
            bt_kwargs['loj'] = False
        else:
            bt_kwargs['wao'] = False
            bt_kwargs['loj'] = True
    else:
        if overlap:
            bt_kwargs['wo'] = True

    with tsv(left) as a, tsv(right) as b:
        out = bedtools.intersect(a=a.name, b=b.name, **bt_kwargs)
        
    bed1_extra_out = bed1_extra.iloc[out[3]].reset_index(drop=True)
    
    if v:        
        out_final = pd.concat([out, bed1_extra_out], axis=1)
    else:
        if outer_join:
            out[4] = out[4].where(out[4] != '.')
            out[7] = out[7].where(out[7] != '.', -1).astype(int)
            
            bed2_extra_out = pd.DataFrame.from_items([
                (name, pd.Series(data=None, index=out.index, dtype=series.dtype))
                for name, series in bed2_extra.iteritems()])
            mask = (out[7] != -1)
            bed2_extra_out.loc[mask, :] = bed2_extra.iloc[out[7][mask]].values
        else:
            bed2_extra_out = bed2_extra.iloc[out[7]].reset_index(drop=True)
        out_final = pd.concat([out, bed1_extra_out, bed2_extra_out], axis=1)
    
    
    outcols = [c + suffixes[0] for c in ['chrom', 'start', 'end', 'index']]
    if not v:
        outcols += [c + suffixes[1] for c in ['chrom', 'start', 'end', 'index']]
    
    if overlap and not v:
        outcols += ['overlap']    
    
    outcols += [c + suffixes[0] for c in bed1_extra_out.columns]
    if not v:
        outcols += [c + suffixes[1] for c in bed2_extra_out.columns]
    
    out_final.columns = outcols
    return out_final
