# -*- coding: utf-8 -*-
import collections

import numpy as np
import pandas as pd

from . import arrops
from ._region import parse_region


def _get_default_colnames():
    return "chrom", "start", "end"


def select(df, region, cols=None):
    """
    Return all genomic intervals in a dataframe that overlap 
    a genomic region.  

    Parameters
    ----------
    df : pandas.DataFrame

    region : UCSC str
        The genomic region to select from the dataframe.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    df : pandas.DataFrame

    """

    ck, sk, ek = _get_default_colnames() if cols is None else cols
    chrom, start, end = parse_region(region)
    if chrom is None:
        raise ValueError("no chromosome detected, check region input")
    if (start is not None) and (end is not None):
        inds = (
            (df.chrom.values == chrom)
            & (df.start.values < end)
            & (df.end.values > start)
        )
    else:
        inds = df.chrom.values == chrom
    return df.iloc[np.where(inds)[0]]


def expand(df, pad, limits=None, side="both", limits_region_col=None, cols=None):
    """
    Expand each interval by a given amount.

    Parameters
    ----------
    df : pandas.DataFrame

    pad : int
        The amount by which the intervals are expanded *on each side*.

    limits : {str: int} or {str: (int, int)}
        The limits of interval expansion. If a single number X if provided,
        the expanded intervals are trimmed to fit into (0, X); if a tuple 
        of numbers is provided (X,Y), the new intervals are trimmed to (X, Y).

    side : str
        Which side to expand, possible values are "left", "right" and "both".

    region_col : str
        The column to select the expansion limits for each interval.
        If None, then use the chromosome column.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    df : pandas.DataFrame
    """

    ck, sk, ek = _get_default_colnames() if cols is None else cols
    limits_region_col = ck if limits_region_col is None else limits_region_col

    if limits:
        lower_limits = {}
        upper_limits = {}
        for k, v in dict(limits).items():
            if isinstance(v, (tuple, list, np.ndarray)):
                lower_limits[k] = v[0]
                upper_limits[k] = v[1]
            elif np.isscalar(v):
                upper_limits[k] = v
                lower_limits[k] = 0
            else:
                raise ValueError("Unknown limit type: {type(v)}")

    if side == "both" or side == "left":
        if limits:
            df[sk] = np.maximum(
                df[limits_region_col].apply(lower_limits.__getitem__, 0),
                df[sk].values - pad,
            )
        else:
            df[sk] = df[sk].values - pad

    if side == "both" or side == "right":
        if limits:
            df[ek] = np.minimum(
                df[limits_region_col].apply(
                    upper_limits.__getitem__, np.iinfo(np.int64).max
                ),
                df[ek] + pad,
            )
        else:
            df[ek] = df[ek] + pad

    return df


def _overlap_intidxs(df1, df2, how="left", keep_order=False, cols1=None, cols2=None):
    """
    Find pairs of overlapping genomic intervals and return the integer
    indices of the overlapping intervals.
    
    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
        
    cols1, cols2 : (str, str, str) or None
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
    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2

    # Switch to integer indices.
    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    # Find overlapping intervals per chromosome.
    df1_groups = df1.groupby(ck1).groups
    df2_groups = df2.groupby(ck2).groups

    overlap_intidxs = []
    for group_keys, df1_group_idxs in df1_groups.items():

        if group_keys not in df2_groups:
            if how == "outer" or how == "left":
                overlap_intidxs_sub = [
                    [
                        df1_group_idxs[:, None],
                        -1 * np.ones_like(df1_group_idxs)[:, None],
                    ]
                ]
                overlap_intidxs.append(np.block(overlap_intidxs_sub))
            continue

        df1_group_idxs = df1_group_idxs.values
        df2_group_idxs = df2_groups[group_keys].values

        df1_group_starts = df1[sk1].values[df1_group_idxs]
        df1_group_ends = df1[ek1].values[df1_group_idxs]
        df2_group_starts = df2[sk2].values[df2_group_idxs]
        df2_group_ends = df2[ek2].values[df2_group_idxs]

        overlap_idxs_loc = arrops.overlap_intervals(
            df1_group_starts, df1_group_ends, df2_group_starts, df2_group_ends,
        )

        # Convert local per-chromosome indices into the
        # indices of the original table.
        overlap_intidxs_sub = [
            [
                df1_group_idxs[overlap_idxs_loc[:, 0]][:, None],
                df2_group_idxs[overlap_idxs_loc[:, 1]][:, None],
            ]
        ]

        if how == "outer" or how == "left":
            no_overlap_ids1 = np.where(
                np.bincount(overlap_idxs_loc[:, 0], minlength=len(df1_group_idxs)) == 0
            )[0]
            overlap_intidxs_sub += [
                [
                    df1_group_idxs[no_overlap_ids1][:, None],
                    -1 * np.ones_like(no_overlap_ids1)[:, None],
                ]
            ]

        if how == "outer" or how == "right":
            no_overlap_ids2 = np.where(
                np.bincount(overlap_idxs_loc[:, 1], minlength=len(df2_group_idxs)) == 0
            )[0]
            overlap_intidxs_sub += [
                [
                    -1 * np.ones_like(no_overlap_ids2)[:, None],
                    df2_group_idxs[no_overlap_ids2][:, None],
                ]
            ]
        overlap_intidxs.append(np.block(overlap_intidxs_sub))

    if how == "outer" or how == "right":
        for group_keys, df2_group_idxs in df2_groups.items():
            if group_keys not in df1_groups:
                overlap_intidxs_sub = [
                    [
                        -1 * np.ones_like(df2_group_idxs)[:, None],
                        df2_group_idxs[:, None],
                    ]
                ]
                overlap_intidxs.append(np.block(overlap_intidxs_sub))

    if len(overlap_intidxs) == 0:
        return np.ndarray(shape=(0, 2), dtype=np.int)
    overlap_intidxs = np.vstack(overlap_intidxs)

    if keep_order:
        order = np.lexsort([overlap_intidxs[:, 1], overlap_intidxs[:, 0]])
        overlap_intidxs = overlap_intidxs[order]

    return overlap_intidxs


def overlap(
    df1,
    df2,
    how="left",
    return_input=True,
    return_index=False,
    return_overlap=False,
    suffixes=["_1", "_2"],
    keep_order=False,
    cols1=None,
    cols2=None,
):

    """
    Find pairs of overlapping genomic intervals.
    
    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
    
    how : {'left', 'right', 'outer', 'inner'}, default 'left'
    
    return

    suffixes : (str, str)
        The suffixes for the columns of the two overlapped sets.
    
    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default 
        values are 'chrom', 'start', 'end'.
    
    Returns
    -------
    df_overlap : pandas.DataFrame
        
    
    """

    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2

    overlap_df_idxs = _overlap_intidxs(
        df1, df2, how=how, cols1=cols1, cols2=cols2, keep_order=keep_order
    )

    # Generate output tables.
    df_index_1 = None
    df_index_2 = None
    if return_index:
        index_col = return_index if isinstance(return_index, str) else "index"
        df_index_1 = pd.DataFrame(
            {index_col + suffixes[0]: df1.index[overlap_df_idxs[:, 0]]}
        )
        df_index_2 = pd.DataFrame(
            {index_col + suffixes[1]: df2.index[overlap_df_idxs[:, 1]]}
        )

    df_overlap = None
    if return_overlap:
        overlap_col = return_overlap if isinstance(return_overlap, str) else "overlap"
        overlap_start = np.maximum(
            df1[sk1].values[overlap_df_idxs[:, 0]],
            df2[sk2].values[overlap_df_idxs[:, 1]],
        )

        overlap_end = np.minimum(
            df1[ek1].values[overlap_df_idxs[:, 0]],
            df2[ek2].values[overlap_df_idxs[:, 1]],
        )

        df_overlap = pd.DataFrame(
            {
                overlap_col + "_" + sk1: overlap_start,
                overlap_col + "_" + ek1: overlap_end,
            }
        )

    df_input_1 = None
    df_input_2 = None
    if return_input == True or str(return_input) == "1" or return_input == "left":
        df_input_1 = df1.iloc[overlap_df_idxs[:, 0]].reset_index(drop=True)
        df_input_1.columns = [c + suffixes[0] for c in df_input_1.columns]
    if return_input == True or str(return_input) == "2" or return_input == "right":
        df_input_2 = df2.iloc[overlap_df_idxs[:, 1]].reset_index(drop=True)
        df_input_2.columns = [c + suffixes[1] for c in df_input_2.columns]

    # Masking non-overlapping regions if using non-inner joins.
    if how != "inner":
        if df_input_1 is not None:
            df_input_1[overlap_df_idxs[:, 0] == -1] = pd.NA
        if df_input_2 is not None:
            df_input_2[overlap_df_idxs[:, 1] == -1] = pd.NA
        if df_index_1 is not None:
            df_index_1[overlap_df_idxs[:, 0] == -1] = pd.NA
        if df_index_2 is not None:
            df_index_2[overlap_df_idxs[:, 1] == -1] = pd.NA
        if df_overlap is not None:
            df_overlap[
                (overlap_df_idxs[:, 0] == -1) | (overlap_df_idxs[:, 1] == -1)
            ] = pd.NA

    out_df = pd.concat(
        [df_index_1, df_input_1, df_index_2, df_input_2, df_overlap], axis="columns"
    )

    return out_df


def cluster(
    df, min_dist=0, out=["input", "cluster"], return_cluster_df=False, cols=None
):
    """
    Cluster overlapping intervals.

    Parameters
    ----------
    df : pandas.DataFrame
    
    min_dist : float or None
        If provided, cluster intervals separated by this distance or less. 
        If None, do not cluster non-overlapping intervals. Using 
        min_dist=0 and min_dist=None will bring different results. 
        bioframe uses semi-open intervals, so interval pairs [0,1) and [1,2)
        do not overlap, but are separated by a distance of 0. Adjacent intervals 
        are not clustered when min_dist=None, but are clustered when min_dist=0.
    
    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.
    
    out : list of str or dict
        A list of requested outputs.
        Can be provided as a dict of {output:column_name} 
        Allowed values: ['input', 'cluster', 'cluster_start', 'cluster_end']
        
    return_cluster_df : bool
        If True, return

    Returns
    -------
    df_clustered : pd.DataFrame

    df_clusters : pd.DataFrame, optional
        A pandas dataframe with coordinates of merged clusters.

    """

    if min_dist is not None:
        if min_dist < 0:
            raise ValueError("min_dist>=0 currently required")

    # Allow users to specify the names of columns containing the interval coordinates.
    ck, sk, ek = _get_default_colnames() if cols is None else cols

    # Switch to integer indices.
    df_index = df.index
    df = df.reset_index(drop=True)

    # Find overlapping intervals per chromosome.
    df_groups = df.groupby(ck).groups

    cluster_ids = np.full(df.shape[0], -1)
    clusters = []
    max_cluster_id = -1

    for group_keys, df_group_idxs in df_groups.items():
        if df_group_idxs.empty:
            continue
        df_group = df.loc[df_group_idxs]

        (
            cluster_ids_group,
            cluster_starts_group,
            cluster_ends_group,
        ) = arrops.merge_intervals(
            df_group[sk].values, df_group[ek].values, min_dist=min_dist
        )

        interval_counts = np.bincount(cluster_ids_group)

        cluster_ids_group += max_cluster_id + 1

        n_clusters = cluster_starts_group.shape[0]
        max_cluster_id += n_clusters

        cluster_ids[df_group_idxs.values] = cluster_ids_group

        ## Storing chromosome names causes a 2x slowdown. :(
        chrom = group_keys
        clusters_group = {
            ck: pd.Series(data=np.full(n_clusters, chrom), dtype=df[ck].dtype),
            sk: cluster_starts_group,
            ek: cluster_ends_group,
            "n_intervals": interval_counts,
        }
        clusters_group = pd.DataFrame(clusters_group)

        clusters.append(clusters_group)

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

    if return_cluster_df:
        return out_df, clusters
    else:
        return out_df


def merge(df, min_dist=0, cols=None):
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
    
    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.
            
    Returns
    -------
    df_merged : pandas.DataFrame
        A pandas dataframe with coordinates of merged clusters.
    """

    if min_dist is not None:
        if min_dist < 0:
            raise ValueError("min_dist>=0 currently required")

    # Allow users to specify the names of columns containing the interval coordinates.
    ck, sk, ek = _get_default_colnames() if cols is None else cols

    # Find overlapping intervals per chromosome.
    df_groups = df.groupby(ck).groups

    clusters = []

    for group_keys, df_group_idxs in df_groups.items():
        if df_group_idxs.empty:
            continue
        df_group = df.loc[df_group_idxs]

        (
            cluster_ids_group,
            cluster_starts_group,
            cluster_ends_group,
        ) = arrops.merge_intervals(
            df_group[sk].values, df_group[ek].values, min_dist=min_dist
        )

        interval_counts = np.bincount(cluster_ids_group)
        n_clusters = cluster_starts_group.shape[0]

        ## Storing chromosome names causes a 2x slowdown. :(
        chrom = group_keys
        clusters_group = {
            ck: pd.Series(data=np.full(n_clusters, chrom), dtype=df[ck].dtype),
            sk: cluster_starts_group,
            ek: cluster_ends_group,
            "n_intervals": interval_counts,
        }
        clusters_group = pd.DataFrame(clusters_group)

        clusters.append(clusters_group)

    clusters = pd.concat(clusters).reset_index(drop=True)

    return clusters


def complement(df, chromsizes=None, cols=None):
    """
    Find genomic regions that are not covered by any interval.

    Parameters
    ----------
    df : pandas.DataFrame
    
    chromsizes : dict

    cols : (str, str, str)
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.
    
    Returns
    -------
    df_complement : numpy.ndarray
    
    """

    # Allow users to specify the names of columns containing the interval coordinates.
    ck, sk, ek = _get_default_colnames() if cols is None else cols

    chromsizes = {} if chromsizes is None else chromsizes

    # Find overlapping intervals per chromosome.
    df_groups = df.groupby(ck).groups

    complements = []

    for group_keys, df_group_idxs in df_groups.items():
        if df_group_idxs.empty:
            continue
        df_group = df.loc[df_group_idxs]

        chrom = group_keys
        if chrom in chromsizes:
            chromsize = chromsizes[chrom]

            if chromsize < np.max(df_group[ek].values):
                raise ValueError("one or more intervals exceed provided chromsize")
            (
                complement_starts_group,
                complement_ends_group,
            ) = arrops.complement_intervals(
                df_group[sk].values, df_group[ek].values, bounds=(0, chromsize),
            )
        else:
            (
                complement_starts_group,
                complement_ends_group,
            ) = arrops.complement_intervals(df_group[sk].values, df_group[ek].values)

        ## Storing chromosome names causes a 2x slowdown. :(
        chrom = group_keys
        complement_group = {
            ck: pd.Series(
                data=np.full(complement_starts_group.shape[0], chrom),
                dtype=df[ck].dtype,
            ),
            sk: complement_starts_group,
            ek: complement_ends_group,
        }
        complement_group = pd.DataFrame(complement_group)

        complements.append(complement_group)

    complements = pd.concat(complements).reset_index(drop=True)

    return complements


def coverage(df1, df2, out=["input", "coverage"], cols1=None, cols2=None):
    """
    Quantify the coverage of intervals from set 1 by intervals from set2. For every interval
     in set 1 find the number of base pairs covered by intervals in set 2. 

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
            
    out : list of str or dict
        A list of requested outputs.
        Can be provided as a dict of {output:column_name} 
        Allowed values: ['input', 'index', 'coverage'].

    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default 
        values are 'chrom', 'start', 'end'.
    
    Returns
    -------
    df_coverage : pandas.DataFrame
    
    """

    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2

    df2_merged = merge(df2, cols=cols2)

    overlap_idxs = overlap(
        df1,
        df2_merged,
        return_index=True,
        return_overlap=True,
        cols1=cols1,
        cols2=cols2,
    )

    overlap_idxs["overlap"] = (
        overlap_idxs["overlap_end"] - overlap_idxs["overlap_start"]
    )

    coverage_sparse_df = overlap_idxs.groupby("index_1").agg({"overlap": "sum"})

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

    out_df = pd.DataFrame(out_df)

    if "input" in out:
        out_df = pd.concat([df1, out_df], axis="columns")

    return out_df


def _closest_intidxs(
    df1,
    df2=None,
    k=1,
    ignore_overlaps=False,
    ignore_upstream=False,
    ignore_downstream=False,
    tie_breaking_col=None,
    cols1=None,
    cols2=None,
):
    """
    For every interval in set 1 find k closest genomic intervals in set2 and
    return their integer indices.
    
    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
        If df2 is None or same object as df1, find closest intervals within the same set.

    k_closest : int
        The number of closest intervals to report.
        
    cols1, cols2 : (str, str, str)
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
    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2

    self_closest = False
    if (df2 is None) or (df2 is df1):
        df2 = df1
        self_closest = True

    # Switch to integer indices.
    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    # Find overlapping intervals per chromosome.
    df1_groups = df1.groupby(ck1).groups
    df2_groups = df2.groupby(ck2).groups
    closest_intidxs = []
    for group_keys, df1_group_idxs in df1_groups.items():
        if group_keys not in df2_groups:
            continue

        df2_group_idxs = df2_groups[group_keys]

        df1_group = df1.loc[df1_group_idxs]
        df2_group = df2.loc[df2_group_idxs]

        tie_arr = None
        if isinstance(tie_breaking_col, str):
            tie_arr = df2_group[tie_breaking_col].values
        elif callable(tie_breaking_col):
            tie_arr = tie_breaking_col(df2_group).values
        else:
            ValueError(
                "tie_breaking_col must be either a column label or "
                "f(DataFrame) -> Series"
            )

        closest_idxs_group = arrops.closest_intervals(
            df1_group[sk1].values,
            df1_group[ek1].values,
            None if self_closest else df2_group[sk2].values,
            None if self_closest else df2_group[ek2].values,
            k=k,
            tie_arr=tie_arr,
            ignore_overlaps=ignore_overlaps,
            ignore_upstream=ignore_upstream,
            ignore_downstream=ignore_downstream,
        )

        # Convert local per-chromosome indices into the
        # indices of the original table.
        closest_idxs_group = np.vstack(
            [
                df1_group_idxs.values[closest_idxs_group[:, 0]],
                df2_group_idxs.values[closest_idxs_group[:, 1]],
            ]
        ).T

        closest_intidxs.append(closest_idxs_group)

    if len(closest_intidxs) == 0:
        return np.ndarray(shape=(0, 2), dtype=np.int)
    closest_intidxs = np.vstack(closest_intidxs)

    return closest_intidxs


def closest(
    df1,
    df2=None,
    k=1,
    ignore_overlaps=False,
    ignore_upstream=False,
    ignore_downstream=False,
    tie_breaking_col=None,
    out=["input", "distance"],
    suffixes=["_1", "_2"],
    cols1=None,
    cols2=None,
):

    """
    For every interval in set 1 find k closest genomic intervals in set2.
    
    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
        If df2 is None, find closest non-identical intervals within the same set.
        
    k : int
        The number of closest intervals to report.
    
    out : list of str or dict
        A list of requested outputs.
        Can be provided as a dict of {output:column_name} 
        Allowed values: ['input', 'index', 'distance', 'have_overlap', 
        'overlap_start', 'overlap_end'].

    suffixes : (str, str)
        The suffixes for the columns of the two sets.
    
    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.
    
    Returns
    -------
    df_closest : pandas.DataFrame
        If no intervals found, returns none.
    
    """
    if k < 1:
        raise ValueError("k>=1 required")

    if df2 is df1:
        raise ValueError(
            "pass df2=None to find closest non-identical intervals within the same set."
        )
    # If finding closest within the same set, df2 now has to be set
    # to df1, so that the rest of the logic works.
    if df2 is None:
        df2 = df1

    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2

    closest_df_idxs = _closest_intidxs(
        df1,
        df2,
        k=k,
        ignore_overlaps=ignore_overlaps,
        ignore_upstream=ignore_upstream,
        ignore_downstream=ignore_downstream,
        tie_breaking_col=tie_breaking_col,
        cols1=cols1,
        cols2=cols2,
    )

    if len(closest_df_idxs) == 0:
        return  # case of no closest intervals

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
            out_df[out["have_overlap"]] = have_overlap
        if "overlap_start" in out:
            out_df[out["overlap_start"]] = np.where(have_overlap, overlap_start, -1)
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


def subtract(df1, df2, cols1=None, cols2=None):

    """
    Generate a new set of genomic intervals by subtracting the second set of genomic intervals from the first.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
    
    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default 
        values are 'chrom', 'start', 'end'.
    
    Returns
    -------
    df_subtracted : pandas.DataFrame
    
    """

    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2

    name_updates = {"chrom_1": "chrom", "overlap_start": "start", "overlap_end": "end"}
    extra_columns_1 = list(np.setdiff1d(df1.columns, [ck1, sk1, ek1]))  # +'_1')
    for i in extra_columns_1:
        name_updates[i + "_1"] = i

    ### loop over chromosomes, then either return the same or subtracted intervals.
    df1_groups = df1.groupby(ck1).groups
    df2_groups = df2.groupby(ck2).groups
    df_subtracted = []
    for group_keys, df1_group_idxs in df1_groups.items():
        df1_group = df1.loc[df1_group_idxs]

        # if nothing to subtract, add original intervals
        if group_keys not in df2_groups:
            df_subtracted.append(df1_group)
            continue

        df2_group_idxs = df2_groups[group_keys]
        df2_group = df2.loc[df2_group_idxs]
        df_subtracted_group = overlap(
            df1_group, complement(df2_group), how="inner", return_overlap=True
        )[list(name_updates)]
        df_subtracted.append(df_subtracted_group.rename(columns=name_updates))
    df_subtracted = pd.concat(df_subtracted)
    return df_subtracted
