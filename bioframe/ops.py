# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import collections

import numpy as np
import pandas as pd

from . import arrops
from .region import parse_region


def bedbisect(bedf, region):
    """Return the span of a block of rows corresponding to
    the genomic region.
    Rows must be sorted by `start` and `end`;
    `chrom` must be grouped, but does not have to be sorted.

    """
    chrom, start, end = parse_region(region)

    lo, hi = arrops._find_block_span(bedf.chrom.values, chrom)

    lo += bedf["end"].values[lo:hi].searchsorted(start, side="right")
    if end is not None:
        hi = lo + bedf["start"].values[lo:hi].searchsorted(end, side="left")
    else:
        hi = None
    return lo, hi


def bedslice(bedf, region):
    """Return a block of rows corresponding to the genomic region.
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
    out = bg2[
        (bg2["chrom1"] == chrom1)
        & (bg2["start1"] >= start1)
        & (bg2["end1"] < end1)
        & (bg2["chrom2"] == chrom2)
        & (bg2["start2"] >= start2)
        & (bg2["end2"] < end2)
    ]
    return out


def expand(df, pad_bp, chromsizes={}, side="both", inplace=False, cols=None):
    ck, sk, ek = ['chrom', 'start', 'end'] if cols is None else cols
    if not inplace:
        df = df.copy()

    if side == "both" or side == "left":
        df[sk] = np.maximum(0, df[sk].values - pad_bp)

    if side == "both" or side == "right":
        if chromsizes:
            df[ek] = np.minimum(
                df[ck].apply(chromsizes.__getitem__, np.iinfo(np.int64).max),
                df[ek] + pad_bp)
        else:
            df[ek] = df[ek] + pad_bp

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
    chroms = kwargs.pop("chroms", None)
    # parallel = kwargs.pop("parallel", False)
    ret_chrom = kwargs.pop("ret_chrom", False)
    map_impl = kwargs.pop("map", map)

    first = tables[0]
    chrom_field = kwargs.pop("chrom_field", first.columns[0])
    if chroms is None:
        chroms = arrops.natsorted(first[chrom_field].unique())

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


def chromsorted(df, by=None, ignore_index=True, chromosomes=None, **kwargs):
    """
    Sort bed-like dataframe by chromosome label in "natural" alphanumeric
    order, followed by any columns specified in ``by``.

    """
    chrom_col = df["chrom"]
    is_categorical = pd.api.types.is_categorical(chrom_col)

    if chromosomes is None:
        if not (is_categorical and chrom_col.cat.ordered):
            dtype = pd.CategoricalDtype(
                arrops.natsorted(chrom_col.unique()), ordered=True
            )
            chrom_col = chrom_col.astype(dtype)
    else:
        dtype = pd.CategoricalDtype(chromosomes, ordered=True)
        chrom_col = chrom_col.astype(dtype)
        missing = df["chrom"].loc[chrom_col.isnull()].unique().tolist()
        if len(missing):
            raise ValueError("Unknown ordering for {}.".format(missing))

    sort_cols = ["chrom"]
    if by is not None:
        if not isinstance(by, list):
            by = [by]
        sort_cols.append(by)

    out = (
        df.assign(chrom=chrom_col)
        .sort_values(sort_cols, **kwargs)
        .reset_index(drop=True)
    )

    if not is_categorical:
        out["chrom"] = out["chrom"].astype(str)

    return out


def make_chromarms(chromsizes, mids, binsize=None, suffixes=("p", "q")):
    """
    Split chromosomes into chromosome arms

    Parameters
    ----------
    chromsizes : pandas.Series
        Series mapping chromosomes to lengths in bp.
    mids : dict-like
        Mapping of chromosomes to midpoint locations.
    binsize : int, optional
        Round midpoints to nearest bin edge for compatibility with a given
        bin grid.
    suffixes : tuple, optional
        Suffixes to name chromosome arms. Defaults to p and q.

    Returns
    -------
    4-column BED-like DataFrame (chrom, start, end, name).
    Arm names are chromosome names + suffix.
    Any chromosome not included in ``mids`` will be omitted.

    """
    chromosomes = [chrom for chrom in chromsizes.index if chrom in mids]

    p_arms = [[chrom, 0, mids[chrom], chrom + suffixes[0]] for chrom in chromosomes]
    if binsize is not None:
        for x in p_arms:
            x[2] = int(round(x[2] / binsize)) * binsize

    q_arms = [
        [chrom, mids[chrom], chromsizes[chrom], chrom + suffixes[1]]
        for chrom in chromosomes
    ]
    if binsize is not None:
        for x in q_arms:
            x[1] = int(round(x[1] / binsize)) * binsize

    interleaved = [*sum(zip(p_arms, q_arms), ())]

    return pd.DataFrame(interleaved, columns=["chrom", "start", "end", "name"])


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
        binedges = np.arange(0, (n_bins + 1)) * binsize
        binedges[-1] = clen
        return pd.DataFrame(
            {"chrom": [chrom] * n_bins, "start": binedges[:-1], "end": binedges[1:]},
            columns=["chrom", "start", "end"],
        )

    bintable = pd.concat(map(_each, chromsizes.keys()), axis=0, ignore_index=True)

    if rel_ids:
        bintable["rel_id"] = bintable.groupby("chrom").cumcount()

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
        raise ValueError("Unknown enzyme name: {}".format(enzyme))

    def _each(chrom):
        seq = bioseq.Seq(str(fasta_records[chrom]))
        cuts = np.r_[0, np.array(cut_finder(seq)) + 1, len(seq)].astype(int)
        n_frags = len(cuts) - 1

        frags = pd.DataFrame(
            {"chrom": [chrom] * n_frags, "start": cuts[:-1], "end": cuts[1:]},
            columns=["chrom", "start", "end"],
        )
        return frags

    return pd.concat(map(_each, chroms), axis=0, ignore_index=True)


def frac_mapped(bintable, fasta_records):
    def _each(bin):
        s = str(fasta_records[bin.chrom][bin.start : bin.end])
        nbases = len(s)
        n = s.count("N")
        n += s.count("n")
        return (nbases - n) / nbases

    return bintable.apply(_each, axis=1)


def frac_gc(bintable, fasta_records, mapped_only=True):
    def _each(chrom_group):
        chrom = chrom_group.name
        seq = fasta_records[chrom]
        gc = []
        for _, bin in chrom_group.iterrows():
            s = str(seq[bin.start : bin.end])
            g = s.count("G")
            g += s.count("g")
            c = s.count("C")
            c += s.count("c")
            nbases = len(s)
            if mapped_only:
                n = s.count("N")
                n += s.count("n")
                nbases -= n
            gc.append((g + c) / nbases if nbases > 0 else np.nan)
        return gc

    out = bintable.groupby("chrom", sort=False).apply(_each)
    return pd.Series(data=np.concatenate(out), index=bintable.index)


def frac_gene_coverage(bintable, mrna):
    if isinstance(mrna, str):
        from .resources import fetch_ucsc_mrna

        mrna = fetch_ucsc_mrna(mrna).rename(
            columns={"tName": "chrom", "tStart": "start", "tEnd": "end"}
        )

    bintable = coverage(
        bintable, 
        mrna, 
        out = {'input':'input', 'count':'gene_count', 'coverage':'gene_coverage'})
    bintable['gene_coverage'] = bintable['gene_coverage'] / (bintable['end'] - bintable['start'])
    
    return bintable


def _overlap_intidxs(df1, df2, **kwargs):
    """
    Find pairs of overlapping genomic intervals and return the integer 
    indices of the overlapping intervals.
    
    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
        
    cols1, cols2 : [str, str, str]
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default 
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    overlap_ids : numpy.ndarray
        The indices of the overlapping genomic intervals in the original 
        dataframes. The 1st column contains the indices of intervals 
        from the 1st set, the 2nd column - the indicies from the 2nd set.
    """

    # Allow users to specify the names of columns containing the interval coordinates.
    ck1, sk1, ek1 = kwargs.get("cols1", ["chrom", "start", "end"])
    ck2, sk2, ek2 = kwargs.get("cols2", ["chrom", "start", "end"])

    # Switch to integer indices.
    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    # Find overlapping intervals per chromosome.
    df1_gb = df1.groupby(ck1)
    df2_gb = df2.groupby(ck2)
    overlap_intidxs = []
    for chrom in df1_gb.groups:
        if chrom not in df2_gb.groups:
            continue

        df1_sub = df1_gb.get_group(chrom)
        df2_sub = df2_gb.get_group(chrom)

        overlap_idxs_loc = arrops.overlap_intervals(
            df1_sub[sk1].values,
            df1_sub[ek1].values,
            df2_sub[sk2].values,
            df2_sub[ek2].values,
        )

        # Convert local per-chromosome indices into the
        # indices of the original table.
        overlap_intidxs_sub = np.vstack(
            [
                df1_gb.groups[chrom][overlap_idxs_loc[:, 0]].values,
                df2_gb.groups[chrom][overlap_idxs_loc[:, 1]].values,
            ]
        ).T

        overlap_intidxs.append(overlap_intidxs_sub)

    overlap_intidxs = np.vstack(overlap_intidxs)

    return overlap_intidxs


def overlap(
    df1,
    df2,
    out=["input", "overlap_start", "overlap_end"],
    suffixes=["_1", "_2"],
    **kwargs
):

    """
    Find pairs of overlapping genomic intervals and return a combined DataFrame.
    
    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
    
    out : list of str or dict
        A list of requested outputs.
        Can be provided as a dict of {output:column_name} 
        Allowed values: ['input', 'index', 'overlap_start', 'overlap_end'].

    suffixes : [str, str]
        The suffixes for the columns of the two overlapped sets.
    
    cols1, cols2 : [str, str, str]
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default 
        values are 'chrom', 'start', 'end'.
    
    Returns
    -------
    df_overlap : pandas.DataFrame
        
    
    """

    ck1, sk1, ek1 = kwargs.get("cols1", ["chrom", "start", "end"])
    ck2, sk2, ek2 = kwargs.get("cols2", ["chrom", "start", "end"])

    overlap_df_idxs = _overlap_intidxs(df1, df2, **kwargs)

    if not isinstance(out, collections.abc.Mapping):
        out = {col: col for col in out}

    out_df = {}
    if "index" in out:
        out_df[out["index"] + suffixes[0]] = df1.index[overlap_df_idxs[:, 0]]
        out_df[out["index"] + suffixes[1]] = df2.index[overlap_df_idxs[:, 1]]
    if "overlap_start" in out:
        out_df[out["overlap_start"]] = np.amax(
            np.vstack(
                [
                    df1[sk1].values[overlap_df_idxs[:, 0]],
                    df2[sk2].values[overlap_df_idxs[:, 1]],
                ]
            ),
            axis=0,
        )
    if "overlap_end" in out:
        out_df[out["overlap_end"]] = np.amin(
            np.vstack(
                [
                    df1[ek1].values[overlap_df_idxs[:, 0]],
                    df2[ek2].values[overlap_df_idxs[:, 1]],
                ]
            ),
            axis=0,
        )
    out_df = pd.DataFrame(out_df)

    if "input" in out:
        df_left = df1.iloc[overlap_df_idxs[:, 0]].reset_index(drop=True)
        df_left.columns = [c + suffixes[0] for c in df_left.columns]
        df_right = df2.iloc[overlap_df_idxs[:, 1]].reset_index(drop=True)
        df_right.columns = [c + suffixes[1] for c in df_right.columns]

        out_df = pd.concat([df_left, df_right, out_df], axis="columns")

    return out_df


def merge(df, min_dist=0, out=["input", "cluster"], **kwargs):
    """
    Merge overlapping intervals.

    Parameters
    ----------
    df : pandas.DataFrame
    
    min_dist : float or None
        If provided, merge intervals separated by this distance or less. 
        If None, do not merge non-overlapping intervals. Using 
        min_dist=0 and min_dist=None will bring different results. 
        bioframe uses semi-open intervals, so interval pairs [0,1) and [1,2)
        do not overlap, but are separated by a distance of 0. Adjacent intervals 
        are not merged when min_dist=None, but are merged when min_dist=0.
    
    cols : [str, str, str]
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.
    
        out : list of str or dict
        A list of requested outputs.
        Can be provided as a dict of {output:column_name} 
        Allowed values: ['input', 'cluster', 'cluster_start', 'cluster_end']
        
    Returns
    -------
    df : numpy.ndarray
    
    clusters : pandas.DataFrame
        A pandas dataframe with coordinates of merged clusters.
    """

    # Allow users to specify the names of columns containing the interval coordinates.
    ck, sk, ek = kwargs.get("cols", ["chrom", "start", "end"])

    # Switch to integer indices.
    df_index = df.index
    df = df.reset_index(drop=True)

    # Find overlapping intervals per chromosome.
    df_gb = df.groupby(ck)

    cluster_ids = np.full(df.shape[0], -1)
    clusters = []
    max_cluster_id = -1

    for chrom, df_chunk in df_gb:
        if df_chunk.empty:
            continue

        (
            cluster_ids_chunk,
            cluster_starts_chunk,
            cluster_ends_chunk,
        ) = arrops.merge_intervals(
            df_chunk[sk].values, df_chunk[ek].values, min_dist=min_dist
        )

        interval_counts = np.bincount(cluster_ids_chunk)

        cluster_ids_chunk += max_cluster_id + 1

        n_clusters = cluster_starts_chunk.shape[0]
        max_cluster_id += n_clusters

        cluster_ids[df_gb.groups[chrom].values] = cluster_ids_chunk

        ## Storing chromosome names causes a 2x slowdown. :(
        clusters_chunk = {
            ck: pd.Series(data=np.full(n_clusters, chrom), dtype=df[ck].dtype),
            sk: cluster_starts_chunk,
            ek: cluster_ends_chunk,
            "count": interval_counts,
        }
        clusters_chunk = pd.DataFrame(clusters_chunk)

        clusters.append(clusters_chunk)

    assert np.all(cluster_ids >= 0)
    clusters = pd.concat(clusters).reset_index(drop=True)

    if not isinstance(out, collections.abc.Mapping):
        out = {col: col for col in out}

    out_df = {}
    if "cluster" in out:
        out_df[out["cluster"]] = cluster_ids
    if "cluster_start" in out:
        out_df[out["cluster_start"]] = clusters[sk].values[cluster_ids]
    if "cluster_end" in out:
        out_df[out["cluster_end"]] = clusters[ek].values[cluster_ids]

    out_df = pd.DataFrame(out_df)

    if "input" in out:
        out_df = pd.concat([df, out_df], axis="columns")

    out_df.set_index(df_index)

    return out_df, clusters


def complement(df, chromsizes={}, **kwargs):
    """
    Find genomic regions that are not covered by at least one of the provided intervals. 

    Parameters
    ----------
    df : pandas.DataFrame
    
    cols : [str, str, str]
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.
    
    Returns
    -------
    df : numpy.ndarray
    
    clusters : pandas.DataFrame
        A pandas dataframe with coordinates of merged clusters.
    """

    # Allow users to specify the names of columns containing the interval coordinates.
    ck, sk, ek = kwargs.get("cols", ["chrom", "start", "end"])

    # Find overlapping intervals per chromosome.
    df_gb = df.groupby(ck)

    complements = []

    for chrom, df_chunk in df_gb:
        if df_chunk.empty:
            continue

        if chrom in chromsizes:
            chromsize = chromsizes[chrom]
            (
                complement_starts_chunk,
                complement_ends_chunk,
            ) = arrops.complement_intervals(
                df_chunk[sk].values, df_chunk[ek].values, bounds=(0, chromsize),
            )
        else:
            (
                complement_starts_chunk,
                complement_ends_chunk,
            ) = arrops.complement_intervals(df_chunk[sk].values, df_chunk[ek].values)

        ## Storing chromosome names causes a 2x slowdown. :(
        complement_chunk = {
            ck: pd.Series(
                data=np.full(complement_starts_chunk.shape[0], chrom),
                dtype=df[ck].dtype,
            ),
            sk: complement_starts_chunk,
            ek: complement_ends_chunk,
        }
        complement_chunk = pd.DataFrame(complement_chunk)

        complements.append(complement_chunk)

    complements = pd.concat(complements).reset_index(drop=True)

    return complements


def coverage(df1, df2, out=["input", "coverage", "count"], **kwargs):
    """
    For every interval in set 1 find the number of overlapping intervals from set 2 and 
    the number of base pairs covered by at least one genomic interval.
    
    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
            
    out : list of str or dict
        A list of requested outputs.
        Can be provided as a dict of {output:column_name} 
        Allowed values: ['input', 'index', 'coverage', 'count'].

    cols1, cols2 : [str, str, str]
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default 
        values are 'chrom', 'start', 'end'.
    
    Returns
    -------
    df_coverage : pandas.DataFrame
    
    """

    ck1, sk1, ek1 = kwargs.get("cols1", ["chrom", "start", "end"])
    ck2, sk2, ek2 = kwargs.get("cols2", ["chrom", "start", "end"])

    df2_merged = merge(df2, **kwargs)

    overlap_idxs = overlap(
        df1,
        df2_merged,
        out=["index", "overlap_start", "overlap_end"],
        cols=[ck2, sk2, ek2],
    )

    overlap_idxs["overlap"] = (
        overlap_idxs["overlap_end"] - overlap_idxs["overlap_start"]
    )

    coverage_sparse_df = (
        overlap_idxs.groupby("index_1")
        .agg({"overlap": "sum", "index_2": "count"})
        .rename(columns={"index_2": "count"})
    )

    # Make an output DataFrame.
    if not isinstance(out, collections.abc.Mapping):
        out = {col: col for col in out}

    out_df = {}

    if "index" in out:
        out_df[out["index"]] = df1.index

    if "coverage" in out:
        out_df[out["coverage"]] = (
            pd.Series(np.zeros_like(df1[sk1]), index=df1.index)
            .add(coverage_sparse_df["overlap"], fill_value=0)
            .astype(df1[sk1].dtype)
        )

    if "count" in out:
        out_df[out["count"]] = (
            pd.Series(np.zeros(df1.shape[0], dtype=np.int64), index=df1.index)
            .add(coverage_sparse_df["count"], fill_value=0)
            .astype(np.int64)
        )

    out_df = pd.DataFrame(out_df)

    if "input" in out:
        out_df = pd.concat([df1, out_df], axis="columns")

    return out_df


def _closest_intidxs(
    df1,
    df2,
    k=1,
    ignore_overlaps=False,
    ignore_upstream=False,
    ignore_downstream=False,
    tie_breaking_col=None,
    **kwargs
):
    """
    For every interval in set 1 find k closest genomic intervals in set2 and
    return their integer indices.
    
    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
        
    k_closest : int
        The number of closest intervals to report.
        
    cols1, cols2 : [str, str, str]
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default 
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    closest_ids : numpy.ndarray
        The indices of the overlapping genomic intervals in the original 
        dataframes. The 1st column contains the indices of intervals 
        from the 1st set, the 2nd column - the indicies from the 2nd set.
    """

    # Allow users to specify the names of columns containing the interval coordinates.
    ck1, sk1, ek1 = kwargs.get("cols1", ["chrom", "start", "end"])
    ck2, sk2, ek2 = kwargs.get("cols2", ["chrom", "start", "end"])

    # Switch to integer indices.
    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    # Find overlapping intervals per chromosome.
    df1_gb = df1.groupby(ck1)
    df2_gb = df2.groupby(ck2)
    closest_intidxs = []
    for chrom in df1_gb.groups:
        if chrom not in df2_gb.groups:
            continue

        df1_chunk = df1_gb.get_group(chrom)
        df2_chunk = df2_gb.get_group(chrom)

        tie_arr = None
        if isinstance(tie_breaking_col, str):
            tie_arr = df2_chunk[tie_breaking_col].values
        elif callable(tie_breaking_col):
            tie_arr = tie_breaking_col(df2_chunk).values
        else:
            ValueError(
                "tie_breaking_col must be either a column label or "
                "f(DataFrame) -> Series"
            )

        closest_idxs_chunk = arrops.closest_intervals(
            df1_chunk[sk1].values,
            df1_chunk[ek1].values,
            df2_chunk[sk2].values,
            df2_chunk[ek2].values,
            k=k,
            tie_arr=tie_arr,
            ignore_overlaps=ignore_overlaps,
            ignore_upstream=ignore_upstream,
            ignore_downstream=ignore_downstream,
        )

        # Convert local per-chromosome indices into the
        # indices of the original table.
        closest_idxs_chunk = np.vstack(
            [
                df1_gb.groups[chrom][closest_idxs_chunk[:, 0]].values,
                df2_gb.groups[chrom][closest_idxs_chunk[:, 1]].values,
            ]
        ).T

        closest_intidxs.append(closest_idxs_chunk)

    closest_intidxs = np.vstack(closest_intidxs)

    return closest_intidxs


def closest(
    df1,
    df2,
    k=1,
    ignore_overlaps=False,
    ignore_upstream=False,
    ignore_downstream=False,
    tie_breaking_col=None,
    out=["input", "distance"],
    suffixes=["_1", "_2"],
    **kwargs
):

    """
    For every interval in set 1 find k closest genomic intervals in set2 and
    return a combined DataFrame.
    
    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
        
    k : int
        The number of closest intervals to report.
    
    out : list of str or dict
        A list of requested outputs.
        Can be provided as a dict of {output:column_name} 
        Allowed values: ['input', 'index', 'distance', 'have_overlap', 
        'overlap_start', 'overlap_end'].

    suffixes : [str, str]
        The suffixes for the columns of the two sets.
    
    cols1, cols2 : [str, str, str]
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.
    
    Returns
    -------
    df_closest : pandas.DataFrame
    
    """

    ck1, sk1, ek1 = kwargs.get("cols1", ["chrom", "start", "end"])
    ck2, sk2, ek2 = kwargs.get("cols2", ["chrom", "start", "end"])

    closest_df_idxs = _closest_intidxs(
        df1,
        df2,
        k=k,
        ignore_overlaps=ignore_overlaps,
        ignore_upstream=ignore_upstream,
        ignore_downstream=ignore_downstream,
        tie_breaking_col=tie_breaking_col,
        **kwargs
    )

    # Make an output DataFrame.
    if not isinstance(out, collections.abc.Mapping):
        out = {col: col for col in out}

    out_df = {}
    if "index" in out:
        out_df[out["index"] + suffixes[0]] = df1.index[closest_df_idxs[:, 0]]
        out_df[out["index"] + suffixes[1]] = df2.index[closest_df_idxs[:, 1]]

    if any([k in out for k in ["have_overlap", "overlap_start", "overlap_end"]]):
        overlap_start = np.amax(
            np.vstack(
                [
                    df1[sk1].values[closest_df_idxs[:, 0]],
                    df2[sk2].values[closest_df_idxs[:, 1]],
                ]
            ),
            axis=0,
        )
        overlap_end = np.amin(
            np.vstack(
                [
                    df1[ek1].values[closest_df_idxs[:, 0]],
                    df2[ek2].values[closest_df_idxs[:, 1]],
                ]
            ),
            axis=0,
        )
        have_overlap = overlap_start < overlap_end

        if "have_overlap" in out:
            out_df[out["overlap_start"]] = have_overlap
        if "overlap_start" in out:
            out_df[out["overlap_start"]] = np.where(
                have_overlap, overlap_start, -1
            )
        if "overlap_end" in out:
            out_df[out["overlap_end"]] = np.where(have_overlap, overlap_end, -1)

    if "distance" in out:
        distance_left = np.maximum(
            0,
            df1[sk1].values[closest_df_idxs[:, 0]]
            - df2[ek2].values[closest_df_idxs[:, 1]],
        )
        distance_right = np.maximum(
            0,
            df2[sk2].values[closest_df_idxs[:, 1]]
            - df1[ek1].values[closest_df_idxs[:, 0]],
        )
        distance = np.amax(np.vstack([distance_left, distance_right]), axis=0)
        out_df[out["distance"]] = distance

    out_df = pd.DataFrame(out_df)

    if "input" in out:
        df_left = df1.iloc[closest_df_idxs[:, 0]].reset_index(drop=True)
        df_left.columns = [c + suffixes[0] for c in df_left.columns]
        df_right = df2.iloc[closest_df_idxs[:, 1]].reset_index(drop=True)
        df_right.columns = [c + suffixes[1] for c in df_right.columns]

        out_df = pd.concat([df_left, df_right, out_df], axis="columns")

    return out_df
