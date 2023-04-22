import numpy as np
import pandas as pd

from .core.specs import _get_default_colnames, _verify_columns
from .core.stringops import parse_region

from .core import arrops
from .core import specs
from .core import construction
from .core import checks

__all__ = [
    "select",
    "select_mask",
    "select_indices",
    "select_labels",
    "expand",
    "overlap",
    "cluster",
    "merge",
    "coverage",
    "closest",
    "subtract",
    "setdiff",
    "count_overlaps",
    "trim",
    "complement",
    "sort_bedframe",
    "assign_view",
]


def select_mask(df, region, cols=None):
    """
    Return boolean mask for all genomic intervals that overlap a query range.

    Parameters
    ----------
    df : pandas.DataFrame

    region : str or tuple
        The genomic region to select from the dataframe in UCSC-style genomic
        region string, or triple (chrom, start, end).

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.

    Returns
    -------
    Boolean array of shape (len(df),)
    """
    ck, sk, ek = _get_default_colnames() if cols is None else cols
    _verify_columns(df, [ck, sk, ek])

    chrom, start, end = parse_region(region)

    if chrom is None:
        raise ValueError("no chromosome detected, check region input")

    if start is None:
        mask = df[ck] == chrom
    else:
        if end is None:
            end = np.inf
        mask = (df[ck] == chrom) & (
            ((df[sk] < end) & (df[ek] > start)) | 
            ((df[sk] == df[ek]) & (df[sk] == start))  # include points at query start
        )
    return mask.to_numpy()


def select_indices(df, region, cols=None):
    """
    Return integer indices of all genomic intervals that overlap a query range.

    Parameters
    ----------
    df : pandas.DataFrame

    region : str or tuple
        The genomic region to select from the dataframe in UCSC-style genomic
        region string, or triple (chrom, start, end).

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.

    Returns
    -------
    1D array of int
    """
    return np.nonzero(select_mask(df, region, cols))[0]


def select_labels(df, region, cols=None):
    """
    Return pandas Index labels of all genomic intervals that overlap a query
    range.

    Parameters
    ----------
    df : pandas.DataFrame

    region : str or tuple
        The genomic region to select from the dataframe in UCSC-style genomic
        region string, or triple (chrom, start, end).

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.

    Returns
    -------
    pandas.Index
    """
    return df.index[select_mask(df, region, cols)]


def select(df, region, cols=None):
    """
    Return all genomic intervals in a dataframe that overlap a genomic region.

    Parameters
    ----------
    df : pandas.DataFrame

    region : str or tuple
        The genomic region to select from the dataframe in UCSC-style genomic
        region string, or triple (chrom, start, end).

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.

    Returns
    -------
    df : pandas.DataFrame

    Notes
    -----
    See :func:`.core.stringops.parse_region()` for more information on region
    formatting.

    See also
    --------
    :func:`select_mask`
    :func:`select_indices`
    :func:`select_labels`
    """
    return df.loc[select_mask(df, region, cols)]


def expand(df, pad=None, scale=None, side="both", cols=None):
    """
    Expand each interval by an amount specified with `pad`.

    Negative values for pad shrink the interval, up to the midpoint.
    Multiplicative rescaling of intervals enabled with scale. Only one of pad
    or scale can be provided. Often followed by :func:`trim()`.

    Parameters
    ----------
    df : pandas.DataFrame

    pad : int, optional
        The amount by which the intervals are additively expanded *on each side*.
        Negative values for pad shrink intervals, but not beyond the interval midpoint.
        Either `pad` or `scale` must be supplied.

    scale : float, optional
        The factor by which to scale intervals multiplicatively on each side, e.g
        ``scale=2`` doubles each interval, ``scale=0`` returns midpoints, and
        ``scale=1`` returns original intervals. Default False.
        Either `pad` or `scale` must be supplied.

    side : str, optional
        Which side to expand, possible values are 'left', 'right' and 'both'.
        Default 'both'.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. Default values are 'chrom', 'start', 'end'.

    Returns
    -------
    df_expanded : pandas.DataFrame

    Notes
    -----
    See :func:`bioframe.trim` for trimming interals after expansion.

    """

    ck, sk, ek = _get_default_colnames() if cols is None else cols
    checks.is_bedframe(df, raise_errors=True, cols=[ck, sk, ek])

    if scale is not None:
        if scale < 0:
            raise ValueError("multiplicative scale must be >=0")
        pads = 0.5 * (scale - 1) * (df[ek].values - df[sk].values)
        types = df.dtypes[[sk, ek]]
    elif pad is not None:
        if not isinstance(pad, int):
            raise ValueError("additive pad must be integer")
        pads = pad
    else:
        raise ValueError("either pad or scale must be supplied")

    df_expanded = df.copy()
    if side == "both" or side == "left":
        df_expanded[sk] = df[sk].values - pads
    if side == "both" or side == "right":
        df_expanded[ek] = df[ek] + pads

    if pad is not None:
        if pad < 0:
            mids = df[sk].values + (0.5 * (df[ek].values - df[sk].values)).astype(int)
            df_expanded[sk] = np.minimum(df_expanded[sk].values, mids)
            df_expanded[ek] = np.maximum(df_expanded[ek].values, mids)
    if scale is not None:
        df_expanded[[sk, ek]] = df_expanded[[sk, ek]].round()
        df_expanded[[sk, ek]] = df_expanded[[sk, ek]].astype(types)

    return df_expanded


def _overlap_intidxs(df1, df2, how="left", cols1=None, cols2=None, on=None):
    """
    Find pairs of overlapping genomic intervals and return the integer
    indices of the overlapping intervals.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.

    how : {'left', 'right', 'outer', 'inner'}, default 'left'
        How to handle the overlaps on the two dataframes.
        left: use the set of intervals in df1
        right: use the set of intervals in df2
        outer: use the union of the set of intervals from df1 and df2
        inner: use intersection of the set of intervals from df1 and df2

    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    on : list or None
        Additional shared columns to consider as separate groups.

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

    # Calculate groups, determined by chrom and on.
    group_list1 = [ck1]
    group_list2 = [ck2]
    if on is not None:
        group_list1 += on
        group_list2 += on
    df1_groups = df1.groupby(group_list1, observed=True, dropna=False).indices

    df2_groups = df2.groupby(group_list2, observed=True, dropna=False).indices
    all_groups = set.union(set(df1_groups), set(df2_groups))

    # Find overlapping intervals per group (determined by chrom and on).
    overlap_intidxs = []
    for group_keys in all_groups:
        df1_group_idxs = (
            df1_groups[group_keys] if (group_keys in df1_groups) else np.array([])
        )
        df2_group_idxs = (
            df2_groups[group_keys] if (group_keys in df2_groups) else np.array([])
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
                [
                    no_overlap_ids1,
                    -1 * np.ones_like(no_overlap_ids1),
                ]
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
                [
                    -1 * np.ones_like(no_overlap_ids2),
                    no_overlap_ids2,
                ]
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
        return np.ndarray(shape=(0, 2), dtype=int)
    overlap_intidxs = np.vstack(overlap_intidxs)

    return overlap_intidxs


def overlap(
    df1,
    df2,
    how="left",
    return_input=True,
    return_index=False,
    return_overlap=False,
    suffixes=("", "_"),
    keep_order=None,
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
        How to handle the overlaps on the two dataframes.
        left: use the set of intervals in df1
        right: use the set of intervals in df2
        outer: use the union of the set of intervals from df1 and df2
        inner: use intersection of the set of intervals from df1 and df2

    return_input : bool
        If True, return columns from input dfs. Default True.

    return_index : bool
        If True, return indicies of overlapping pairs as two new columns
        ('index'+suffixes[0] and 'index'+suffixes[1]). Default False.

    return_overlap : bool
        If True, return overlapping intervals for the overlapping pairs
        as two additional columns (`overlap_start`, `overlap_end`). Default False.

    suffixes : (str, str)
        The suffixes for the columns of the two overlapped sets.

    keep_order : bool, optional
        If True and how='left', sort the output dataframe to preserve the order
        of the intervals in df1. Cannot be used with how='right'/'outer'/'inner'.
        Default True for how='left', and None otherwise.

    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    on : list or None
        List of additional shared columns to consider as separate groups
        when considering overlaps. A common use would be passing on=['strand'].
        Default is None.

    Returns
    -------
    df_overlap : pandas.DataFrame

    """

    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2
    checks.is_bedframe(df1, raise_errors=True, cols=[ck1, sk1, ek1])
    checks.is_bedframe(df2, raise_errors=True, cols=[ck2, sk2, ek2])

    if (how == "left") and (keep_order is None):
        keep_order = True
    if (how != "left") and (keep_order is True):
        raise ValueError("keep_order=True only allowed for how='left'")

    if on is None:
        on_list = []
    else:
        if not isinstance(on, list):
            raise ValueError("on=[] must be None or list")
        if (ck1 in on) or (ck2 in on):
            raise ValueError("on=[] should not contain chromosome colnames")
        _verify_columns(df1, on)
        _verify_columns(df2, on)
        on_list = on

    overlap_df_idxs = _overlap_intidxs(
        df1,
        df2,
        how=how,
        cols1=cols1,
        cols2=cols2,
        on=on,
    )

    # Generate output tables.
    index_col = return_index if isinstance(return_index, str) else "index"
    index_col_1 = index_col + suffixes[0]
    index_col_2 = index_col + suffixes[1]
    df_index_1 = pd.DataFrame({index_col_1: df1.index[overlap_df_idxs[:, 0]]})
    df_index_2 = pd.DataFrame({index_col_2: df2.index[overlap_df_idxs[:, 1]]})

    df_overlap = None
    if return_overlap:
        overlap_col = return_overlap if isinstance(return_overlap, str) else "overlap"
        overlap_col_sk1 = overlap_col + "_" + sk1
        overlap_col_ek1 = overlap_col + "_" + ek1

        overlap_start = np.maximum(
            df1[sk1].values[overlap_df_idxs[:, 0]],
            df2[sk2].values[overlap_df_idxs[:, 1]],
        )

        overlap_end = np.minimum(
            df1[ek1].values[overlap_df_idxs[:, 0]],
            df2[ek2].values[overlap_df_idxs[:, 1]],
        )

        df_overlap = pd.DataFrame(
            {overlap_col_sk1: overlap_start, overlap_col_ek1: overlap_end}
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
        df_index_1[overlap_df_idxs[:, 0] == -1] = None
        df_index_1 = df_index_1.astype({index_col_1: pd.Int64Dtype()})
        df_index_2[overlap_df_idxs[:, 1] == -1] = None
        df_index_2 = df_index_2.astype({index_col_2: pd.Int64Dtype()})

        if df_input_1 is not None:
            df_input_1[overlap_df_idxs[:, 0] == -1] = None
            df_input_1 = df_input_1.astype(
                {
                    (sk1 + suffixes[0]): pd.Int64Dtype(),
                    (ek1 + suffixes[0]): pd.Int64Dtype(),
                }
            )
        if df_input_2 is not None:
            df_input_2[overlap_df_idxs[:, 1] == -1] = None
            df_input_2 = df_input_2.astype(
                {
                    (sk2 + suffixes[1]): pd.Int64Dtype(),
                    (ek2 + suffixes[1]): pd.Int64Dtype(),
                }
            )
        if df_overlap is not None:
            df_overlap[
                (overlap_df_idxs[:, 0] == -1) | (overlap_df_idxs[:, 1] == -1)
            ] = None
            df_overlap = df_overlap.astype(
                {overlap_col_sk1: pd.Int64Dtype(), (overlap_col_ek1): pd.Int64Dtype()}
            )

    out_df = pd.concat(
        [df_index_1, df_input_1, df_index_2, df_input_2, df_overlap], axis="columns"
    )
    if keep_order:
        out_df = out_df.sort_values([index_col_1, index_col_2])

    if not return_index:
        out_df.drop([index_col_1, index_col_2], axis=1, inplace=True)

    out_df.reset_index(drop=True, inplace=True)
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
    Cluster overlapping intervals into groups.

    Can return numeric ids for these groups (with `return_cluster_ids`=True) and/or
    their genomic coordinates (with `return_cluster_intervals`=True).
    Also see :func:`merge()`, which discards original intervals and returns a new set.

    Parameters
    ----------
    df : pandas.DataFrame

    min_dist : float or None
        If provided, cluster intervals separated by this distance or less.
        If ``None``, do not cluster non-overlapping intervals.
        Since bioframe uses semi-open intervals, interval pairs [0,1) and [1,2)
        do not overlap, but are separated by a distance of 0. Such adjacent intervals
        are not clustered when ``min_dist=None``, but are clustered when ``min_dist=0``.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.

    on : None or list
        List of column names to perform clustering on independently, passed as an argument
        to df.groupby before clustering. Default is ``None``.
        An example useage would be to pass ``on=['strand']``.

    return_input : bool
        If True, return input

    return_cluster_ids : bool
        If True, return ids for clusters

    return_cluster_invervals : bool
        If True, return clustered interval the original interval belongs to

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
        if not isinstance(on, list):
            raise ValueError("on=[] must be None or list")
        if ck in on:
            raise ValueError("on=[] should not contain chromosome colnames")
        _verify_columns(df, on)
        group_list += on
    df_groups = df.groupby(group_list, observed=True).groups

    cluster_ids = np.full(df.shape[0], -1)
    clusters = []
    max_cluster_id = -1

    for group_keys, df_group_idxs in df_groups.items():
        if pd.isna(pd.Series(group_keys)).any():
            continue
        if df_group_idxs.empty:
            continue
        df_group = df.loc[df_group_idxs]
        (
            cluster_ids_group,
            cluster_starts_group,
            cluster_ends_group,
        ) = arrops.merge_intervals(
            df_group[sk].values.astype(int),
            df_group[ek].values.astype(int),
            min_dist=min_dist,
        )

        interval_counts = np.bincount(cluster_ids_group)
        cluster_ids_group += max_cluster_id + 1
        n_clusters = cluster_starts_group.shape[0]
        max_cluster_id += n_clusters

        cluster_ids[df_group_idxs.values] = cluster_ids_group

        ## Storing chromosome names causes a 2x slowdown. :(
        if isinstance(group_keys, str):
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

    df_nans = pd.isnull(df[[sk, ek] + group_list]).any(axis=1)
    if df_nans.sum() > 0:
        cluster_ids[df_nans.values] = (
            max_cluster_id + 1 + np.arange(np.sum(df_nans.values))
        )
        clusters.append(df.loc[df_nans])

    clusters = pd.concat(clusters).reset_index(drop=True)
    if df_nans.sum() > 0:
        clusters = clusters.astype({sk: pd.Int64Dtype(), ek: pd.Int64Dtype()})
    assert np.all(cluster_ids >= 0)

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

    This returns a new dataframe of genomic intervals, which have the genomic coordinates
    of the interval cluster groups from the input dataframe. Also :func:`cluster()`, which
    returns the assignment of intervals to clusters prior to merging.

    Parameters
    ----------
    df : pandas.DataFrame

    min_dist : float or None
        If provided, merge intervals separated by this distance or less.
        If None, do not merge non-overlapping intervals. Using
        ``min_dist=0`` and ``min_dist=None`` will bring different results.
        bioframe uses semi-open intervals, so interval pairs [0,1) and [1,2)
        do not overlap, but are separated by a distance of 0. Adjacent intervals
        are not merged when ``min_dist=None``, but are merged when ``min_dist=0``.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.

    on : None or list
        List of column names to perform clustering on independently, passed as an argument
        to df.groupby before clustering. Default is None.
        An example useage would be to pass ``on=['strand']``.

    Returns
    -------
    df_merged : pandas.DataFrame
        A pandas dataframe with coordinates of merged clusters.

    Notes
    -------
    Resets index.

    """

    if min_dist is not None:
        if min_dist < 0:
            raise ValueError("min_dist>=0 currently required")

    # Allow users to specify the names of columns containing the interval coordinates.
    ck, sk, ek = _get_default_colnames() if cols is None else cols
    checks.is_bedframe(df, raise_errors=True, cols=[ck, sk, ek])

    df = df.copy()
    df.reset_index(inplace=True, drop=True)

    # Find overlapping intervals for groups specified by on=[] (default on=None)
    group_list = [ck]
    if on is not None:
        if not isinstance(on, list):
            raise ValueError("on=[] must be None or list")
        if ck in on:
            raise ValueError("on=[] should not contain chromosome colnames")
        _verify_columns(df, on)
        group_list += on
    df_groups = df.groupby(group_list, observed=True).groups

    clusters = []

    for group_keys, df_group_idxs in df_groups.items():
        if pd.isna(pd.Series(group_keys)).any():
            continue
        if df_group_idxs.empty:
            continue

        df_group = df.loc[df_group_idxs]
        (
            cluster_ids_group,
            cluster_starts_group,
            cluster_ends_group,
        ) = arrops.merge_intervals(
            df_group[sk].values.astype(int),
            df_group[ek].values.astype(int),
            min_dist=min_dist
            # df_group[sk].values, df_group[ek].values, min_dist=min_dist
        )
        interval_counts = np.bincount(cluster_ids_group)
        n_clusters = cluster_starts_group.shape[0]

        ## Storing chromosome names causes a 2x slowdown. :(
        if isinstance(group_keys, str):
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

    df_nans = pd.isnull(df[[sk, ek] + group_list]).any(axis=1)
    df_has_nans = df_nans.sum()
    if df_has_nans:
        nan_intervals = pd.DataFrame(
            [pd.NA] * df_has_nans,
            columns=["n_intervals"],
            index=df.loc[df_nans].index,
        )
        clusters.append(
            pd.concat(
                [df.loc[df_nans], nan_intervals],
                axis=1,
            )
        )

    clusters = pd.concat(clusters).reset_index(drop=True)
    if df_has_nans:
        clusters = clusters.astype(
            {sk: pd.Int64Dtype(), ek: pd.Int64Dtype(), "n_intervals": pd.Int64Dtype()}
        )

    # reorder cluster columns to have chrom,start,end first
    clusters_names = list(clusters.keys())
    clusters = clusters[
        [ck, sk, ek] + [col for col in clusters_names if col not in [ck, sk, ek]]
    ]

    return clusters


def coverage(
    df1,
    df2,
    suffixes=("", "_"),
    return_input=True,
    cols1=None,
    cols2=None,
):
    """
    Quantify the coverage of intervals from 'df1' by intervals from 'df2'.

    For every interval in 'df1' find the number of base pairs covered by intervals in 'df2'.
    Note this only quantifies whether a basepair in 'df1' was covered, as 'df2' is merged
    before calculating coverage.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.

    suffixes : (str, str)
        The suffixes for the columns of the two overlapped sets.

    return_input : bool
        If True, return input as well as computed coverage

    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    df_coverage : pandas.DataFrame

    Notes
    ------
    Resets index.

    """

    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2

    df1.reset_index(inplace=True, drop=True)

    df2_merged = merge(df2, cols=cols2)

    df_overlap = overlap(
        df1,
        df2_merged,
        how="left",
        suffixes=suffixes,
        keep_order=True,
        return_index=True,
        return_overlap=True,
        cols1=cols1,
        cols2=cols2,
    )

    df_overlap["overlap"] = df_overlap["overlap_end"] - df_overlap["overlap_start"]

    out_df = (
        pd.DataFrame(
            df_overlap.groupby("index" + suffixes[0])
            .agg({"overlap": "sum"})["overlap"]
            .astype(df1[sk1].dtype)
        )
        .rename(columns={"overlap": "coverage"})
        .reset_index(drop=True)
    )

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
    direction_col=None,
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
        The second column is filled with -1 for those intervals in the 1st
        set with no closest 2nd set interval.
    """

    # Allow users to specify the names of columns containing the interval coordinates.
    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2

    self_closest = False
    if (df2 is None) or (df2 is df1):
        if len(df1) == 1:
            raise ValueError(
                "df1 must have more than one interval to find closest non-identical interval"
            )
        df2 = df1
        self_closest = True

    # Switch to integer indices.
    df1 = df1.reset_index(drop=True)
    df2 = df2.reset_index(drop=True)

    # Find overlapping intervals per chromosome.
    df1_groups = df1.groupby(ck1, observed=True).groups
    df2_groups = df2.groupby(ck2, observed=True).groups
    closest_intidxs = []
    for group_keys, df1_group_idxs in df1_groups.items():
        if group_keys not in df2_groups:
            #
            closest_idxs_group = np.vstack(
                [
                    df1_group_idxs,
                    -1 * np.ones_like(df1_group_idxs),
                ]
            ).T
            closest_intidxs.append(closest_idxs_group)
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

        # Verify and construct the direction_arr (convert from pandas string column to bool array)
        # TODO: should we add checks that it's valid "strand"?
        direction_arr = None
        if direction_col is None:
            direction_arr = np.ones(len(df1_group), dtype=np.bool_)
        else:
            direction_arr = (df1_group[direction_col].values != "-") # both "+" and "." keep orientation by genomic coordinate

        # Calculate closest intervals with arrops:
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
            direction=direction_arr
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
    direction_col=None,
    tie_breaking_col=None,
    return_input=True,
    return_index=False,
    return_distance=True,
    return_overlap=False,
    suffixes=("", "_"),
    cols1=None,
    cols2=None,
):
    """
    For every interval in dataframe `df1` find k closest genomic intervals in dataframe `df2`.

    Currently, we are not taking the feature strands into account for filtering.
    However, the strand can be used for definition of upstream/downstream of the feature (direction).

    Note that, unless specified otherwise, overlapping intervals are considered as closest.
    When multiple intervals are located at the same distance, the ones with the lowest index
    in `df2` are returned.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.
        If `df2` is None, find closest non-identical intervals within the same set.

    k : int
        The number of the closest intervals to report.

    ignore_overlaps : bool
        If True, ignore overlapping intervals and return the closest non-overlapping interval.

    ignore_upstream : bool
        If True, ignore intervals in `df2` that are upstream of intervals in `df1`,
        relative to the reference strand or the strand specified by direction_col.

    ignore_downstream : bool
        If True, ignore intervals in `df2` that are downstream of intervals in `df1`,
        relative to the reference strand or the strand specified by direction_col.

    direction_col : str
        Name of direction column that will set upstream/downstream orientation for each feature.
        The column should contain bioframe-compliant strand values ("+", "-", ".").

    tie_breaking_col : str
        A column in `df2` to use for breaking ties when multiple intervals
        are located at the same distance. Intervals with *lower* values will
        be selected.

    return_input : bool
        If True, return input

    return_index : bool
        If True, return indices

    return_distance : bool
        If True, return distances. Returns zero for overlaps.

    return_overlap : bool
        If True, return columns: 'have_overlap', 'overlap_start', and 'overlap_end'.
        Fills df_closest['overlap_start'] and df['overlap_end']
        with None if non-overlapping. Default False.

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

    Notes
    -----
    By default, direction is defined by the reference genome: everything with
    smaller coordinate is considered upstream, everything with larger coordinate
    is considered downstream.

    If ``direction_col`` is provided, upstream/downstream are relative to the
    direction column in ``df1``, i.e. features marked "+" and "." strand will
    define upstream and downstream as above, while features marked "-" have
    upstream and downstream reversed: smaller coordinates are downstream and
    larger coordinates are upstream.
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
    checks.is_bedframe(df1, raise_errors=True, cols=[ck1, sk1, ek1])
    checks.is_bedframe(df2, raise_errors=True, cols=[ck2, sk2, ek2])

    closest_df_idxs = _closest_intidxs(
        df1,
        df2,
        k=k,
        ignore_overlaps=ignore_overlaps,
        ignore_upstream=ignore_upstream,
        ignore_downstream=ignore_downstream,
        direction_col=direction_col,
        tie_breaking_col=tie_breaking_col,
        cols1=cols1,
        cols2=cols2,
    )
    na_mask = closest_df_idxs[:, 1] == -1

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
        df_index_2[na_mask] = pd.NA

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

        df_overlap = pd.DataFrame(
            {
                "have_overlap": have_overlap,
                "overlap_start": np.where(have_overlap, overlap_start, None),
                "overlap_end": np.where(have_overlap, overlap_end, None),
            },
        )
        df_overlap = df_overlap.astype(
            {
                "have_overlap": pd.BooleanDtype(),
                "overlap_start": pd.Int64Dtype(),
                "overlap_end": pd.Int64Dtype(),
            }
        )
        df_overlap[na_mask] = pd.NA

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
        df_distance = pd.DataFrame({"distance": distance}, dtype=pd.Int64Dtype())
        df_distance[na_mask] = pd.NA

    df_input_1 = None
    df_input_2 = None
    if return_input is True or str(return_input) == "1" or return_input == "left":
        df_input_1 = df1.iloc[closest_df_idxs[:, 0]].reset_index(drop=True)
        df_input_1.columns = [c + suffixes[0] for c in df_input_1.columns]
    if return_input is True or str(return_input) == "2" or return_input == "right":
        df_input_2 = df2.iloc[closest_df_idxs[:, 1]].reset_index(drop=True)
        df_input_2.columns = [c + suffixes[1] for c in df_input_2.columns]
        df_input_2 = df_input_2.astype(
            {sk2 + suffixes[1]: pd.Int64Dtype(), ek2 + suffixes[1]: pd.Int64Dtype()}
        )
        df_input_2[na_mask] = pd.NA

    out_df = pd.concat(
        [df_index_1, df_input_1, df_index_2, df_input_2, df_overlap, df_distance],
        axis="columns",
    )

    return out_df


def subtract(
    df1,
    df2,
    return_index=False,
    suffixes=("", "_"),
    cols1=None,
    cols2=None,
):
    """
    Generate a new set of genomic intervals by subtracting the second set of genomic intervals from the first.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.

    return_index : bool
        Whether to return the indices of the original intervals ('index'+suffixes[0]),
        and the indices of any sub-intervals split by subtraction ('sub_index'+suffixes[1]).
        Default False.

    suffixes : (str,str)
        Suffixes for returned indices. Only alters output if return_index is True.
        Default ("","_").

    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    df_subtracted : pandas.DataFrame

    Notes
    -----
    Resets index, drops completely subtracted (null) intervals, and casts to pd.Int64Dtype().

    """

    ck1, sk1, ek1 = _get_default_colnames() if cols1 is None else cols1
    ck2, sk2, ek2 = _get_default_colnames() if cols2 is None else cols2

    name_updates = {
        ck1 + suffixes[0]: ck1,
        "overlap_" + sk1: sk1,
        "overlap_" + ek1: ek1,
    }
    extra_columns_1 = [i for i in list(df1.columns) if i not in [ck1, sk1, ek1]]
    for i in extra_columns_1:
        name_updates[i + suffixes[0]] = i
    if return_index:
        name_updates["index" + suffixes[0]] = "index" + suffixes[0]
        name_updates["index" + suffixes[1]] = "complement_index" + suffixes[1]

    all_chroms = np.unique(
        list(pd.unique(df1[ck1].dropna())) + list(pd.unique(df2[ck2].dropna()))
    )
    if len(all_chroms) == 0:
        raise ValueError("No chromosomes remain after dropping nulls")

    df_subtracted = overlap(
        df1,
        complement(
            df2, view_df={i: np.iinfo(np.int64).max for i in all_chroms}, cols=cols2
        ).astype({sk2: pd.Int64Dtype(), ek2: pd.Int64Dtype()}),
        how="left",
        suffixes=suffixes,
        return_index=return_index,
        return_overlap=True,
        keep_order=True,
        cols1=cols1,
        cols2=cols2,
    )[list(name_updates)]
    df_subtracted.rename(columns=name_updates, inplace=True)
    df_subtracted = df_subtracted.iloc[~pd.isna(df_subtracted[sk1].values)]
    df_subtracted.reset_index(drop=True, inplace=True)

    if return_index:
        inds = df_subtracted["index" + suffixes[0]].values
        comp_inds = df_subtracted["complement_index" + suffixes[1]].copy()  # .values
        for i in np.unique(inds):
            comp_inds[inds == i] -= comp_inds[inds == i].min()
        df_subtracted["sub_index" + suffixes[1]] = comp_inds.copy()
        df_subtracted.drop(columns=["complement_index" + suffixes[1]], inplace=True)
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
        genomic intervals, provided separately for each dataframe.
        The default values are 'chrom', 'start', 'end'.

    on : None or list
        Additional column names to perform clustering on independently, passed as an argument
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


def count_overlaps(
    df1,
    df2,
    suffixes=("", "_"),
    return_input=True,
    cols1=None,
    cols2=None,
    on=None,
):
    """
    Count number of overlapping genomic intervals.

    Parameters
    ----------
    df1, df2 : pandas.DataFrame
        Two sets of genomic intervals stored as a DataFrame.

    suffixes : (str, str)
        The suffixes for the columns of the two overlapped sets.

    return_input : bool
        If True, return columns from input dfs. Default True.

    cols1, cols2 : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    on : list
        List of additional shared columns to consider as separate groups
        when considering overlaps. A common use would be passing on=['strand'].
        Default is None.

    Returns
    -------
    df_counts : pandas.DataFrame

    Notes
    -------
    Resets index.

    """

    df1.reset_index(inplace=True, drop=True)

    df_counts = overlap(
        df1,
        df2,
        how="left",
        return_input=False,
        keep_order=True,
        suffixes=suffixes,
        return_index=True,
        on=on,
        cols1=cols1,
        cols2=cols2,
    )

    out_df = pd.DataFrame(
        df_counts.groupby(["index" + suffixes[0]])["index" + suffixes[1]]
        .count()
        .values,
        columns=["count"],
    )

    if return_input:
        out_df = pd.concat([df1, out_df], axis="columns")
    return out_df


def trim(
    df,
    view_df=None,
    df_view_col=None,
    view_name_col="name",
    return_view_columns=False,
    cols=None,
    cols_view=None,
):
    """
    Trim each interval to fall within regions specified in the viewframe 'view_df'.

    Intervals that fall outside of view regions are replaced with nulls.
    If no 'view_df' is provided, intervals are truncated at zero to avoid
        negative values.

    Parameters
    ----------
    df : pandas.DataFrame

    view_df : None or pandas.DataFrame
        View specifying region start and ends for trimming. Attempts to
        convert dictionary and pd.Series formats to viewFrames.

    df_view_col : str or None
        The column of 'df' used to specify view regions.
        The associated region in 'view_df' is then used for trimming.
        If None, :func:'bioframe.ops.assign_view' will be used to assign view regions.
        If no 'view_df' is provided, uses the 'chrom' column, df[cols[0]].
        Default None.

    view_name_col : str
        Column of df with region names. Default 'name'.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.

    cols_view : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals in the view. The default values are 'chrom', 'start', 'end'.

    Returns
    -------
    df_trimmed : pandas.DataFrame

    """

    ck, sk, ek = _get_default_colnames() if cols is None else cols
    _verify_columns(df, [ck, sk, ek])
    df_columns = list(df.columns)
    df_trimmed = df.copy()
    inferred_view = False

    if view_df is None:
        df_view_col = ck
        view_df = {
            i: np.iinfo(np.int64).max
            for i in set(df[df_view_col].copy().dropna().values)
        }
        inferred_view = True

    ckv, skv, ekv = _get_default_colnames() if cols_view is None else cols_view
    view_df = construction.make_viewframe(
        view_df, view_name_col=view_name_col, cols=[ckv, skv, ekv]
    ).rename(columns=dict(zip([ckv, skv, ekv], [ck, sk, ek])))

    if inferred_view:
        view_name_col = ck
    elif df_view_col is None:
        if _verify_columns(df_trimmed, ["view_region"], return_as_bool=True):
            raise ValueError("column view_region already exists in input df")
        df_view_col = "view_region"

        df_trimmed = assign_view(
            df_trimmed,
            view_df,
            drop_unassigned=False,
            df_view_col=df_view_col,
            view_name_col=view_name_col,
            cols=cols,
            cols_view=cols,
        )
    else:
        _verify_columns(df_trimmed, [df_view_col])
        checks.is_cataloged(
            df_trimmed,
            view_df,
            raise_errors=True,
            df_view_col=df_view_col,
            view_name_col=view_name_col,
        )

    df_trimmed = df_trimmed.merge(
        view_df,
        how="left",
        left_on=df_view_col,
        right_on=view_name_col,
        suffixes=("", "_view"),
    )

    unassigned_intervals = pd.isnull(df_trimmed[df_view_col].values)
    if unassigned_intervals.any():
        df_trimmed.loc[unassigned_intervals, [ck, sk, ek]] = pd.NA
        df_trimmed.astype({sk: pd.Int64Dtype(), ek: pd.Int64Dtype()})

    lower_vector = df_trimmed[sk + "_view"].values
    upper_vector = df_trimmed[ek + "_view"].values

    df_trimmed[sk].clip(lower=lower_vector, upper=upper_vector, inplace=True)
    df_trimmed[ek].clip(lower=lower_vector, upper=upper_vector, inplace=True)

    if return_view_columns:
        return df_trimmed
    else:
        return df_trimmed[df_columns]


def complement(df, view_df=None, view_name_col="name", cols=None, cols_view=None):
    """
    Find genomic regions in a viewFrame 'view_df' that are not covered by any interval in the dataFrame 'df'.

    First assigns intervals in 'df' to region in 'view_df', splitting intervals in 'df' as necessary.

    Parameters
    ----------
    df : pandas.DataFrame

    view_df : pandas.Dataframe
        If none, attempts to infer the view from chroms (i.e. df[cols[0]]).

    view_name_col : str
        Name of column in view_df with unique reigon names. Default 'name'.

    cols : (str, str, str)
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.

    cols_view : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals in the view. The default values are 'chrom', 'start', 'end'.

    Returns
    -------
    df_complement : pandas.DataFrame

    Notes
    ------
    Discards null intervals in input, and df_complement has regular int dtype.

    """

    ### TODO add on=, so can do strand-specific complements...

    # Allow users to specify the names of columns containing the interval coordinates.
    ck, sk, ek = _get_default_colnames() if cols is None else cols
    _verify_columns(df, [ck, sk, ek])

    if view_df is None:
        view_df = {i: np.iinfo(np.int64).max for i in set(df[ck].dropna().values)}

    ckv, skv, ekv = _get_default_colnames() if cols_view is None else cols_view
    view_df = construction.make_viewframe(
        view_df, view_name_col=view_name_col, cols=[ckv, skv, ekv]
    ).rename(columns=dict(zip([ckv, skv, ekv], [ck, sk, ek])))

    # associate intervals to regions, required to enable single interval from df to
    # overlap multiple intervals in view_df. note this differs from the goal of assign_view.
    new_intervals = overlap(
        view_df,
        df,
        return_overlap=True,
        how="inner",
        suffixes=("", "_df"),
        cols1=cols,
        cols2=cols,
    )
    new_intervals = new_intervals[
        [ck, "overlap_" + sk, "overlap_" + ek, view_name_col]
    ].copy()
    new_intervals.rename(
        columns={
            "overlap_" + sk: sk,
            "overlap_" + ek: ek,
            view_name_col: "view_region",
        },
        inplace=True,
    )
    df = new_intervals
    checks.is_cataloged(
        df,
        view_df,
        raise_errors=True,
        df_view_col="view_region",
        view_name_col=view_name_col,
    )

    # Find overlapping intervals per region.
    df_groups = df.groupby("view_region").groups
    all_groups = sorted(set(view_df[view_name_col]))

    complements = []
    for group_key in all_groups:
        region_interval = view_df.loc[view_df[view_name_col] == group_key]
        region_chrom, region_start, region_end = region_interval[[ck, sk, ek]].values[0]

        if group_key not in df_groups:
            complement_group = region_interval.copy().rename(
                columns={view_name_col: "view_region"}
            )
            complements.append(pd.DataFrame(complement_group))
            continue

        df_group_idxs = df_groups[group_key].values
        df_group = df.loc[df_group_idxs]

        (complement_starts_group, complement_ends_group,) = arrops.complement_intervals(
            df_group[sk].values.astype(int),
            df_group[ek].values.astype(int),
            bounds=(region_start, region_end),
        )

        # Storing chromosome names causes a 2x slowdown. :(
        complement_group = {
            ck: pd.Series(
                data=np.full(complement_starts_group.shape[0], region_chrom),
                dtype=df[ck].dtype,
            ),
            sk: complement_starts_group,
            ek: complement_ends_group,
            "view_region": group_key,
        }
        complement_group = pd.DataFrame(complement_group)

        complements.append(complement_group)

    complements = pd.concat(complements).reset_index(drop=True)

    return complements


def sort_bedframe(
    df,
    view_df=None,
    reset_index=True,
    df_view_col=None,
    view_name_col="name",
    cols=None,
    cols_view=None,
):
    """
    Sorts a bedframe 'df'.

    If 'view_df' is not provided, sorts by ``cols`` (e.g. "chrom", "start", "end").
    If 'view_df' is provided and 'df_view_col' is not provided, uses
    :func:`bioframe.ops.assign_view` with ``df_view_col='view_region'``
    to assign intervals to the view regions with the largest overlap and then sorts.
    If 'view_df' and 'df_view_col' are both provided, checks if the latter
    are cataloged in 'view_name_col', and then sorts.

    df : pandas.DataFrame
        Valid bedframe.

    view_df : pandas.DataFrame | dict-like
        Valid input to make a viewframe. When it is dict-like
        :func:'bioframe.make_viewframe' will be used to convert
        to viewframe. If view_df is not provided df is sorted by chrom and start.

    reset_index : bool
        Default True.

    df_view_col: None | str
        Column from 'df' used to associate intervals with view regions.
        The associated region in 'view_df' is then used for sorting.
        If None, :func:'bioframe.assign_view' will be used to assign view regions.
        Default None.

    view_name_col: str
        Column from view_df with names of regions.
        Default `name`.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.

    cols_view : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals in the view. The default values are 'chrom', 'start', 'end'.

    Returns
    -------
    out_df : sorted bedframe

    Notes
    -------
        df_view_col is currently returned as an ordered categorical

    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    if not checks.is_bedframe(df, cols=cols):
        raise ValueError("not a valid bedframe, cannot sort")

    out_df = df.copy()
    if view_df is None:
        out_df.sort_values([ck1, sk1, ek1], inplace=True)

    else:
        ckv, skv, ekv = _get_default_colnames() if cols_view is None else cols_view
        view_df = construction.make_viewframe(
            view_df, view_name_col=view_name_col, cols=[ckv, skv, ekv]
        ).rename(columns=dict(zip([ckv, skv, ekv], [ck1, sk1, ek1])))

        if df_view_col is None:
            if _verify_columns(out_df, ["view_region"], return_as_bool=True):
                raise ValueError("column view_region already exists in input df")
            df_view_col = "view_region"
            out_df = assign_view(
                out_df,
                view_df,
                df_view_col=df_view_col,
                view_name_col=view_name_col,
                cols=cols,
                cols_view=cols,
            )

        else:
            if not _verify_columns(out_df, [df_view_col], return_as_bool=True):
                raise ValueError(
                    "column 'df_view_col' not in input df, cannot sort by view"
                )
            if not checks.is_cataloged(
                out_df[pd.isna(out_df[df_view_col].values) == False],
                view_df,
                df_view_col=df_view_col,
                view_name_col=view_name_col,
            ):
                raise ValueError(
                    "intervals in df not cataloged in view_df, cannot sort by view"
                )

        view_cat = pd.CategoricalDtype(
            categories=view_df[view_name_col].values, ordered=True
        )
        out_df[df_view_col] = out_df[df_view_col].astype({df_view_col: view_cat})
        out_df.sort_values([df_view_col, ck1, sk1, ek1], inplace=True)

    # make sure no columns get appended and dtypes are preserved
    out_df = out_df[df.columns].astype(df.dtypes)

    if reset_index:
        out_df.reset_index(inplace=True, drop=True)

    return out_df


def assign_view(
    df,
    view_df,
    drop_unassigned=False,
    df_view_col="view_region",
    view_name_col="name",
    cols=None,
    cols_view=None,
):
    """
    Associates genomic intervals in bedframe ``df`` with regions in viewframe ``view_df``, based on their largest overlap.

    Parameters
    ----------
    df : pandas.DataFrame

    view_df : pandas.DataFrame
        ViewFrame specifying region start and ends for assignment. Attempts to
        convert dictionary and pd.Series formats to viewFrames.

    drop_unassigned : bool
        If True, drop intervals in df that do not overlap a region in the view.
        Default False.

    df_view_col : str
        The column of ``df`` used to specify view regions.
        The associated region in view_df is then used for trimming.
        If no view_df is provided, uses the chrom column, ``df[cols[0]]``.
        Default "view_region".

    view_name_col : str
        Column of ``view_df`` with region names. Default 'name'.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.

    cols_view : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals in the view. The default values are 'chrom', 'start', 'end'.

    Returns
    -------
    out_df : dataframe with an associated view region for each interval in ``out_df[view_name_col]``.

    Notes
    -------
    Resets index.

    """

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    df = df.copy()
    df.reset_index(inplace=True, drop=True)

    checks.is_bedframe(df, raise_errors=True, cols=cols)
    view_df = construction.make_viewframe(
        view_df, view_name_col=view_name_col, cols=cols_view
    )

    overlap_view = overlap(
        df,
        view_df,
        how="left",
        suffixes=("", "_view"),
        return_overlap=True,
        keep_order=False,
        return_index=True,
        cols1=cols,
        cols2=cols_view,
    )

    overlap_columns = overlap_view.columns
    overlap_view["overlap_length"] = (
        overlap_view["overlap_" + ek1] - overlap_view["overlap_" + sk1]
    )

    out_df = (
        overlap_view.sort_values("overlap_length", ascending=False)
        .drop_duplicates("index", keep="first")
        .sort_values("index")
    )

    out_df.rename(columns={view_name_col + "_view": df_view_col}, inplace=True)

    if drop_unassigned:
        out_df = out_df.iloc[pd.isna(out_df).any(axis=1).values == 0, :]
    out_df.reset_index(inplace=True, drop=True)

    return_cols = list(df.columns) + [df_view_col]

    return out_df[return_cols]
