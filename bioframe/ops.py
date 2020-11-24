import numpy as np
import pandas as pd
import collections

from . import arrops
from .region import parse_region, regions_add_name_column

_rc = {
    'colnames':{
        'chrom':'chrom',
        'start':'start',
        'end':'end'
    }
}


def _get_default_colnames():
    return _rc['colnames']['chrom'], _rc['colnames']['start'], _rc['colnames']['end']


class update_default_colnames:
    def __init__(self, new_colnames):
        self._old_colnames = dict(_rc['colnames'])
        if isinstance(new_colnames, collections.Iterable):
            if len(new_colnames) != 3:
                raise ValueError(
                    'Please, specify new columns using a list of '
                    '3 strings or a dict!')
            (_rc['colnames']['chrom'],
            _rc['colnames']['start'],
            _rc['colnames']['end']) = new_colnames
        elif isinstance(new_colnames, collections.Mapping):
            _rc['colnames'].update({k:v for k,v in new_colnames.items()
                                    if k in ['chrom', 'start', 'end']})
        else:
            raise ValueError(
                    'Please, specify new columns using a list of '
                    '3 strings or a dict!')

    def __enter__(self):
        return self

    def __exit__(self, *args):
        _rc['colnames'] = self._old_colnames


def _verify_columns(df, colnames):
    """
    df: pandas.DataFrame

    colnames: list of columns
    """
    if not set(colnames).issubset(df.columns):
        raise ValueError(
            ", ".join(set(colnames).difference(set(df.columns)))
            + " not in keys of df.columns"
        )


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


def _overlap_intidxs(
    df1, df2, how="left", keep_order=False, cols1=None, cols2=None, on=None
):
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

    on : list or None
        Additional shared columns to consider as separate groups

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
    _verify_columns(df1, [ck1, sk1, ek1])
    _verify_columns(df2, [ck2, sk2, ek2])

    # Switch to integer indices.
    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    # Find overlapping intervals per chromosome.
    group_list1 = [ck1]
    group_list2 = [ck2]
    if on is not None:
        if type(on) is not list:
            raise ValueError("on=[] must be None or list")
        if (ck1 in on) or (ck2 in on):
            raise ValueError("on=[] should not contain chromosome colnames")
        _verify_columns(df1, on)
        _verify_columns(df2, on)
        group_list1 += on
        group_list2 += on
    df1_groups = df1.groupby(group_list1).groups
    df2_groups = df2.groupby(group_list2).groups
    all_groups = sorted(
        set.union(set(df1_groups), set(df2_groups))
    )  ### breaks if any of the groupby elements are pd.NA...
    # all_groups = list(set.union(set(df1_groups), set(df2_groups))) ### disagrees with pyranges order so a test fails...

    overlap_intidxs = []
    for group_keys in all_groups:
        df1_group_idxs = (
            df1_groups[group_keys].values
            if (group_keys in df1_groups)
            else np.array([])
        )
        df2_group_idxs = (
            df2_groups[group_keys].values
            if (group_keys in df2_groups)
            else np.array([])
        )
        overlap_intidxs_sub = []

        both_groups_nonempty = (df1_group_idxs.size > 0) and (df2_group_idxs.size > 0)

        if both_groups_nonempty:
            overlap_idxs_loc = arrops.overlap_intervals(
                df1[sk1].values[df1_group_idxs],
                df1[ek1].values[df1_group_idxs],
                df2[sk2].values[df2_group_idxs],
                df2[ek2].values[df2_group_idxs],
            )

            # Convert local per-chromosome indices into the
            # indices of the original table.
            overlap_intidxs_sub += [
                [
                    df1_group_idxs[overlap_idxs_loc[:, 0]],
                    df2_group_idxs[overlap_idxs_loc[:, 1]],
                ]
            ]

        if how in ["outer", "left"] and df1_group_idxs.size > 0:
            if both_groups_nonempty:
                no_overlap_ids1 = df1_group_idxs[
                    np.where(
                        np.bincount(
                            overlap_idxs_loc[:, 0], minlength=len(df1_group_idxs)
                        )
                        == 0
                    )[0]
                ]
            else:
                no_overlap_ids1 = df1_group_idxs

            overlap_intidxs_sub += [
                [no_overlap_ids1, -1 * np.ones_like(no_overlap_ids1),]
            ]

        if how in ["outer", "right"] and df2_group_idxs.size > 0:
            if both_groups_nonempty:
                no_overlap_ids2 = df2_group_idxs[
                    np.where(
                        np.bincount(
                            overlap_idxs_loc[:, 1], minlength=len(df2_group_idxs)
                        )
                        == 0
                    )[0]
                ]
            else:
                no_overlap_ids2 = df2_group_idxs

            overlap_intidxs_sub += [
                [-1 * np.ones_like(no_overlap_ids2), no_overlap_ids2,]
            ]
        if overlap_intidxs_sub:
            overlap_intidxs.append(
                np.block(
                    [
                        [idxs[:, None] for idxs in idxs_pair]
                        for idxs_pair in overlap_intidxs_sub
                    ]
                )
            )

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
    suffixes=("_1", "_2"),
    keep_order=False,
    cols1=None,
    cols2=None,
    on=None,
):

    """
    Find pairs of overlapping genomic intervals.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.

    how : {'left', 'right', 'outer', 'inner'}, default 'left'

    return_input : bool
        If True, return columns from input dfs. Default True.

    return_index : bool
        If True, return indicies of overlapping pairs. Default False.

    return_overlap
        If True, return overlapping intervals for the overlapping pairs. Default False.

    suffixes : (str, str)
        The suffixes for the columns of the two overlapped sets.

    keep_order : bool
        << to be documented >>

    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    on : list
        List of column names to perform clustering on indepdendently, passed as an argument
        to df.groupby when considering overlaps. Default is ['chrom'], which must match the first name
        from cols. Examples for additional columns include 'strand'.

    Returns
    -------
    df_overlap : pandas.DataFrame

    """

    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2

    overlap_df_idxs = _overlap_intidxs(
        df1, df2, how=how, cols1=cols1, cols2=cols2, keep_order=keep_order, on=on,
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
    if return_input is True or str(return_input) == "1" or return_input == "left":
        df_input_1 = df1.iloc[overlap_df_idxs[:, 0]].reset_index(drop=True)
        df_input_1.columns = [c + suffixes[0] for c in df_input_1.columns]
    if return_input is True or str(return_input) == "2" or return_input == "right":
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
    df,
    min_dist=0,
    cols=None,
    on=None,
    return_input=True,
    return_cluster_ids=True,
    return_cluster_intervals=True,
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

    on : None or list
        List of column names to perform clustering on indepdendently, passed as an argument
        to df.groupby before clustering. Default is None. An example use would be on=['strand'].

    return_input : bool
        If True, return input

    return_cluster_ids : bool
        If True, return ids for clusters

    return_cluster_invervals : bool
        If True, return clustered interval the original interval belongs to

    return_cluster_df : bool
        If True, return df_clusters

    Returns
    -------
    df_clustered : pd.DataFrame

    """
    if min_dist is not None:
        if min_dist < 0:
            raise ValueError("min_dist>=0 currently required")
    # Allow users to specify the names of columns containing the interval coordinates.
    ck, sk, ek = _get_default_colnames() if cols is None else cols
    _verify_columns(df, [ck, sk, ek])

    # Switch to integer indices.
    df_index = df.index
    df = df.reset_index(drop=True)

    # Find overlapping intervals for groups specified by ck1 and on=[] (default on=None)
    group_list = [ck]
    if on is not None:
        if type(on) is not list:
            raise ValueError("on=[] must be None or list")
        if ck in on:
            raise ValueError("on=[] should not contain chromosome colnames")
        _verify_columns(df, on)
        group_list += on
    df_groups = df.groupby(group_list).groups

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
        if type(group_keys) is str:
            group_keys = tuple((group_keys,))
        clusters_group = {}
        for col in group_list:
            clusters_group[col] = pd.Series(
                data=np.full(n_clusters, group_keys[group_list.index(col)]),
                dtype=df[col].dtype,
            )
        clusters_group[sk] = cluster_starts_group
        clusters_group[ek] = cluster_ends_group
        clusters_group["n_intervals"] = interval_counts
        clusters_group = pd.DataFrame(clusters_group)
        clusters.append(clusters_group)

    assert np.all(cluster_ids >= 0)
    clusters = pd.concat(clusters).reset_index(drop=True)
    # reorder cluster columns to have chrom,start,end first
    clusters_names = list(clusters.keys())
    clusters = clusters[
        [ck, sk, ek] + [col for col in clusters_names if col not in [ck, sk, ek]]
    ]

    out_df = {}
    if return_cluster_ids:
        out_df["cluster"] = cluster_ids
    if return_cluster_intervals:
        out_df["cluster_start"] = clusters[sk].values[cluster_ids]
        out_df["cluster_end"] = clusters[ek].values[cluster_ids]

    out_df = pd.DataFrame(out_df)

    if return_input:
        out_df = pd.concat([df, out_df], axis="columns")

    out_df.set_index(df_index)

    return out_df


def merge(df, min_dist=0, cols=None, on=None):
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

    on : list
        List of column names to consider separately for merging, passed as an argument
        to df.groupby before merging. Default is ['chrom'], which must match the first name
        from cols. Examples for additional columns include 'strand'.

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
    _verify_columns(df, [ck, sk, ek])

    # Find overlapping intervals for groups specified by on=[] (default on=None)
    group_list = [ck]
    if on is not None:
        if type(on) is not list:
            raise ValueError("on=[] must be None or list")
        if ck in on:
            raise ValueError("on=[] should not contain chromosome colnames")
        _verify_columns(df, on)
        group_list += on
    df_groups = df.groupby(group_list).groups

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
        if type(group_keys) is str:
            group_keys = tuple((group_keys,))
        clusters_group = {}
        for col in group_list:
            clusters_group[col] = pd.Series(
                data=np.full(n_clusters, group_keys[group_list.index(col)]),
                dtype=df[col].dtype,
            )
        clusters_group[sk] = cluster_starts_group
        clusters_group[ek] = cluster_ends_group
        clusters_group["n_intervals"] = interval_counts
        clusters_group = pd.DataFrame(clusters_group)

        clusters.append(clusters_group)

    clusters = pd.concat(clusters).reset_index(drop=True)
    # reorder cluster columns to have chrom,start,end first
    clusters_names = list(clusters.keys())
    clusters = clusters[
        [ck, sk, ek] + [col for col in clusters_names if col not in [ck, sk, ek]]
    ]

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

    infer_chromsizes = (chromsizes is None)

    # Find overlapping intervals per chromosome.
    df_groups = df.groupby(ck).groups

    if infer_chromsizes:
        all_groups = sorted(set(df_groups))
    else:
        if not set(df_groups).issubset(set(chromsizes.keys())):
            raise ValueError(
                'Chromsizes are missing some chromosomes from the input interval table.')
        all_groups = sorted(set(chromsizes.keys()))

    complements = []

    for group_keys in all_groups:
        # this is a stub for potential on argument
        chrom = group_keys

        if group_keys not in df_groups:
            complement_group = {
                ck: pd.Series(
                    data=[chrom],
                    dtype=df[ck].dtype,
                ),
                sk: 0,
                ek: chromsizes[chrom],
            }

            complements.append(pd.DataFrame(complement_group))
            continue

        df_group_idxs = df_groups[group_keys].values
        df_group = df.loc[df_group_idxs]

        if infer_chromsizes:
            chromsize = np.iinfo(np.int64).max
        else:
            chromsize = chromsizes[chrom]

        if chromsize < np.max(df_group[ek].values):
            raise ValueError("one or more intervals exceed provided chromsize")
        (
            complement_starts_group,
            complement_ends_group,
        ) = arrops.complement_intervals(
            df_group[sk].values, df_group[ek].values, bounds=(0, chromsize),
        )

        # Storing chromosome names causes a 2x slowdown. :(
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


def coverage(df1, df2, return_input=True, cols1=None, cols2=None):
    """
    Quantify the coverage of intervals from set 1 by intervals from set2. For every interval
     in set 1 find the number of base pairs covered by intervals in set 2.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.

    return_input : bool
        If True, return input as well as computed coverage

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
        how="left",
        return_index=True,
        return_overlap=True,
        cols1=cols1,
        cols2=cols2,
    )

    overlap_idxs["overlap"] = (
        overlap_idxs["overlap_end"] - overlap_idxs["overlap_start"]
    )

    coverage_sparse_df = overlap_idxs.groupby("index_1").agg({"overlap": "sum"})

    out_df = {}
    out_df["coverage"] = (
        pd.Series(np.zeros_like(df1[sk1]), index=df1.index)
        .add(coverage_sparse_df["overlap"], fill_value=0)
        .astype(df1[sk1].dtype)
    )
    out_df = pd.DataFrame(out_df)

    if return_input:
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
    return_input=True,
    return_index=False,
    return_distance=True,
    return_overlap=False,
    suffixes=("_1", "_2"),
    cols1=None,
    cols2=None,
):

    """
    For every interval in set 1 find k closest genomic intervals in set 2.
    Note that, unless specified otherwise, overlapping intervals are considered
    as closest. When multiple intervals are located at the same distance, the
    ones with the lowest index in df2 are chosen.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
        If df2 is None, find closest non-identical intervals within the same set.

    k : int
        The number of closest intervals to report.

    ignore_overlaps : bool
        If True, return the closest non-overlapping interval.

    ignore_upstream : bool
        If True, ignore intervals in df2 that are upstream of intervals in df1.

    ignore_downstream : bool
        If True, ignore intervals in df2 that are downstream of intervals in df1.

    tie_breaking_col : str
        A column in df2 to use for breaking ties when multiple intervals
        are located at the same distance. Intervals with *lower* values will
        be selected.

    return_input : bool
        If True, return input

    return_index : bool
        If True, return indices

    return_distance : bool
        If True, return distances. Returns zero for overlaps.

    return_overlap = False,
        If True, return columns: have_overlap, overlap_start, and overlap_end.
        Fills df_closest['overlap_start'] and df['overlap_end'] with pd.NA if non-overlapping.

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

    # Generate output tables.
    df_index_1 = None
    df_index_2 = None
    if return_index:
        index_col = return_index if isinstance(return_index, str) else "index"
        df_index_1 = pd.DataFrame(
            {index_col + suffixes[0]: df1.index[closest_df_idxs[:, 0]]}
        )
        df_index_2 = pd.DataFrame(
            {index_col + suffixes[1]: df2.index[closest_df_idxs[:, 1]]}
        )

    df_overlap = None
    if return_overlap:
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

        df_overlap = pd.DataFrame({
            "have_overlap" : have_overlap,
            "overlap_start" : np.where(have_overlap, overlap_start, pd.NA),
            "overlap_end": np.where(have_overlap, overlap_end, pd.NA)
        })

    df_distance = None
    if return_distance:
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
        df_distance = pd.DataFrame({
            "distance" : distance
        })


    df_input_1 = None
    df_input_2 = None
    if return_input is True or str(return_input) == "1" or return_input == "left":
        df_input_1 = df1.iloc[closest_df_idxs[:, 0]].reset_index(drop=True)
        df_input_1.columns = [c + suffixes[0] for c in df_input_1.columns]
    if return_input is True or str(return_input) == "2" or return_input == "right":
        df_input_2 = df2.iloc[closest_df_idxs[:, 1]].reset_index(drop=True)
        df_input_2.columns = [c + suffixes[1] for c in df_input_2.columns]

    out_df = pd.concat([
        df_index_1,
        df_input_1,
        df_index_2,
        df_input_2,
        df_overlap,
        df_distance], axis="columns")

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

    name_updates = {ck1 + "_1": "chrom", "overlap_start": "start", "overlap_end": "end"}
    extra_columns_1 = [i for i in list(df1.columns) if i not in [ck1, sk1, ek1]]
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


def setdiff(df1, df2, cols1=None, cols2=None, on=None):
    """
    Generate a new dataframe of genomic intervals by removing any interval from the
    first dataframe that overlaps an interval from the second dataframe.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as DataFrames.

    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    on : None or list
        Additional column names to perform clustering on indepdendently, passed as an argument
        to df.groupby when considering overlaps and must be present in both dataframes.
        Examples for additional columns include 'strand'.

    Returns
    -------
    df_setdiff : pandas.DataFrame

    """
    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2

    df_overlapped = _overlap_intidxs(
        df1, df2, how="inner", cols1=cols1, cols2=cols2, on=on
    )
    inds_non_overlapped = np.setdiff1d(np.arange(len(df1)), df_overlapped[:, 0])
    df_setdiff = df1.iloc[inds_non_overlapped]
    return df_setdiff


def split(
    df,
    points,
    cols=None,
    cols_points=None,
    add_names=False,
    suffixes=["_left", "_right"],
):
    """
    Generate a new dataframe of genomic intervals by splitting each interval from the
    first dataframe that overlaps an interval from the second dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        Genomic intervals stored as a DataFrame.

    points : pandas.DataFrame or dict
        If pandas.DataFrame, a set of genomic positions specified in columns 'chrom', 'pos'.
        Names of cols can be overwridden by cols_points.
        If dict, mapping of chromosomes to positions.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.


    Returns
    -------
    df_split : pandas.DataFrame

    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    ck2, sk2 = ("chrom", "pos") if cols_points is None else cols_points

    name_updates = {ck1 + "_1": "chrom", "overlap_"+sk1: "start", "overlap_"+ek1: "end"}
    if add_names:
        name_updates["index_2"] = "index_2"
        return_index = True
    else:
        return_index = False
    extra_columns_1 = [i for i in list(df.columns) if i not in [ck1, sk1, ek1]]
    for i in extra_columns_1:
        name_updates[i + "_1"] = i

    if isinstance(points, dict):
        points = pd.DataFrame.from_dict(points, orient="index", columns=[sk2])
        points.reset_index(inplace=True)
        points.rename(columns={"index": "chrom"}, inplace=True)
    elif not isinstance(points, pd.DataFrame):
        raise ValueError("points must be a dict or pd.Dataframe")

    points["start"] = points[sk2]
    points["end"] = points[sk2]
    all_chroms = set(df[ck1].unique()).union(df[ck2].unique()) 
    all_chroms = {c:np.iinfo(np.int64).max for c in all_chroms}
    df_split = overlap(
        df,
        complement(points, chromsizes=all_chroms, cols=(ck2, 'start', 'end')),
        how="inner",
        cols1=cols,
        cols2=(ck2, "start", "end"),
        return_overlap=True,
        return_index=return_index,
    )[list(name_updates)]
    df_split.rename(columns=name_updates, inplace=True)

    if add_names:
        df_split = regions_add_name_column(df_split)
        sides = np.mod(df_split["index_2"].values, 2).astype(int)  # .astype(str)
        df_split["name"] = df_split["name"].values + np.array(suffixes)[sides]
        df_split.drop(columns=["index_2"])

    return df_split


def count_overlaps(
    df1, df2, cols1=None, cols2=None, on=None,
):

    """
    Count number of overlapping genomic intervals.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.

    how : {'left', 'right', 'outer', 'inner'}, default 'left'

    return_input : bool
        If True, return columns from input dfs. Default True.

    suffixes : (str, str)
        The suffixes for the columns of the two overlapped sets.

    keep_order : bool
        << to be documented >>

    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    on : list
        List of column names to check overlap on indepdendently, passed as an argument
        to df.groupby when considering overlaps. Default is None. Examples for additional columns include 'strand'.

    Returns
    -------
    df_counts : pandas.DataFrame

    """

    df_counts = overlap(
        df1,
        df2,
        how="left",
        return_input=False,
        keep_order=True,
        return_index=True,
        on=on,
        cols1=cols1,
        cols2=cols2,
    )
    df_counts = pd.concat(
        [
            df1,
            pd.DataFrame(
                df_counts.groupby(["index_1"])["index_2"].count().values,
                columns=["count"],
            ),
        ],
        axis=1,
        names=["count"],
    )

    return df_counts
