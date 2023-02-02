import re
import warnings

import numpy as np
import pandas as pd


def natsort_key(s, _NS_REGEX=re.compile(r"(\d+)", re.U)):
    return tuple([int(x) if x.isdigit() else x for x in _NS_REGEX.split(s) if x])


def natsorted(iterable):
    return sorted(iterable, key=natsort_key)


def argnatsort(array):
    array = np.asarray(array)
    if not len(array):
        return np.array([], dtype=int)
    cols = tuple(zip(*(natsort_key(x) for x in array)))
    return np.lexsort(cols[::-1])  # numpy's lexsort is ass-backwards


def _find_block_span(arr, val):
    """Find the first and the last occurence + 1 of the value in the array."""
    # it can be done via bisection, but for now BRUTE FORCE
    block_idxs = np.where(arr == val)[0]
    lo, hi = block_idxs[0], block_idxs[-1] + 1
    return lo, hi


def interweave(a, b):
    """
    Interweave two arrays.

    Parameters
    ----------
    a, b : numpy.ndarray
        Arrays to interweave, must have the same length/

    Returns
    -------
    out : numpy.ndarray
        Array of interweaved values from a and b.

    Notes
    -----
    From https://stackoverflow.com/questions/5347065/interweaving-two-numpy-arrays
    """
    out = np.empty((a.size + b.size,), dtype=a.dtype)
    out[0::2] = a
    out[1::2] = b
    return out


def sum_slices(arr, starts, ends):
    """
    Calculate sums of slices of an array.

    Parameters
    ----------
    arr : numpy.ndarray
    starts : numpy.ndarray
        Starts for each slice
    ends : numpy.ndarray
        Stops for each slice

    Returns
    -------
    sums : numpy.ndarray
        Sums of the slices.
    """
    sums = np.add.reduceat(arr, interweave(starts, ends))[::2]
    sums[starts == ends] = 0
    return sums


def arange_multi(starts, stops=None, lengths=None):
    """
    Create concatenated ranges of integers for multiple start/length.

    Parameters
    ----------
    starts : numpy.ndarray
        Starts for each range
    stops : numpy.ndarray
        Stops for each range
    lengths : numpy.ndarray
        Lengths for each range. Either stops or lengths must be provided.

    Returns
    -------
    concat_ranges : numpy.ndarray
        Concatenated ranges.

    Notes
    -----
    See the following illustrative example:

    starts = np.array([1, 3, 4, 6])
    stops = np.array([1, 5, 7, 6])

    print arange_multi(starts, lengths)
    >>> [3 4 4 5 6]

    From: https://codereview.stackexchange.com/questions/83018/vectorized-numpy-version-of-arange-with-multiple-start-stop

    """

    if (stops is None) == (lengths is None):
        raise ValueError("Either stops or lengths must be provided!")

    if lengths is None:
        lengths = stops - starts

    if np.isscalar(starts):
        starts = np.full(len(stops), starts)

    # Repeat start position index length times and concatenate
    cat_start = np.repeat(starts, lengths)

    # Create group counter that resets for each start/length
    cat_counter = np.arange(lengths.sum()) - np.repeat(
        lengths.cumsum() - lengths, lengths
    )

    # Add group counter to group specific starts
    cat_range = cat_start + cat_counter

    return cat_range


def _check_overlap(starts1, ends1, starts2, ends2, closed=False):
    """
    Take pairs of intervals and test if each pair has an overlap.

    Parameters
    ----------
    starts1, ends1, starts2, ends2 : numpy.ndarray
        Interval coordinates. All four arrays must have the same size.
        Warning: if provided as pandas.Series, indices will be ignored.

    closed : bool
        If True then treat intervals as closed and accept single-point overlaps.

    Returns
    -------
    have_overlap : numpy.ndarray
        A boolean array where the i-th element says if the i-th interval in set 1
        overlaps the i-th interval in set 2.
    """

    if not (starts1.size == ends1.size == starts2.size == ends2.size):
        raise ValueError("All four input arrays must have the same size.")

    if closed:
        return (starts1 <= ends2) & (starts2 <= ends1)
    else:
        return (starts1 < ends2) & (starts2 < ends1)


def _size_overlap(starts1, ends1, starts2, ends2):
    """
    Take pairs of intervals and return the length of an overlap in each pair.

    Parameters
    ----------
    starts1, ends1, starts2, ends2 : numpy.ndarray
        Interval coordinates. All four arrays must have the same size.
        Warning: if provided as pandas.Series, indices will be ignored.

    Returns
    -------
    overlap_size : numpy.ndarray
        An array where the i-th element contains the length of an overlap between
        the i-th interval in set 1 and the i-th interval in set 2.
        0 if the intervals overlap by a single point, -1 if they do not overlap.
    """

    overlap_size = np.minimum(ends1, ends2) - np.maximum(starts1, starts2)
    overlap_size[overlap_size < 0] = -1
    return overlap_size


def _overlap_intervals_legacy(starts1, ends1, starts2, ends2, closed=False, sort=False):
    """
    Take two sets of intervals and return the indices of pairs of overlapping intervals.

    Parameters
    ----------
    starts1, ends1, starts2, ends2 : numpy.ndarray
        Interval coordinates. Warning: if provided as pandas.Series, indices
        will be ignored.

    closed : bool
        If True, then treat intervals as closed and report single-point overlaps.

    Returns
    -------
    overlap_ids : numpy.ndarray
        An Nx2 array containing the indices of pairs of overlapping intervals.
        The 1st column contains ids from the 1st set, the 2nd column has ids
        from the 2nd set.

    """

    for vec in [starts1, ends1, starts2, ends2]:
        if issubclass(type(vec), pd.core.series.Series):
            warnings.warn(
                "One of the inputs is provided as pandas.Series and its index "
                "will be ignored.",
                SyntaxWarning,
            )

    starts1 = np.asarray(starts1)
    ends1 = np.asarray(ends1)
    starts2 = np.asarray(starts2)
    ends2 = np.asarray(ends2)

    # Concatenate intervals lists
    n1 = len(starts1)
    n2 = len(starts2)
    starts = np.concatenate([starts1, starts2])
    ends = np.concatenate([ends1, ends2])

    # Encode interval ids as 1-based,
    # negative ids for the 1st set, positive ids for 2nd set
    ids = np.concatenate([-np.arange(1, n1 + 1), np.arange(1, n2 + 1)])

    # Sort all intervals together
    order = np.lexsort([ends, starts])
    starts, ends, ids = starts[order], ends[order], ids[order]

    # Find interval overlaps
    match_starts = np.arange(0, n1 + n2)
    match_ends = np.searchsorted(starts, ends, "right" if closed else "left")

    # Ignore self-overlaps
    match_mask = match_ends > match_starts + 1
    match_starts, match_ends = match_starts[match_mask], match_ends[match_mask]

    # Restore
    overlap_ids = np.vstack(
        [
            np.repeat(ids[match_starts], match_ends - match_starts - 1),
            ids[arange_multi(match_starts + 1, match_ends)],
        ]
    ).T

    # Drop same-set overlaps
    overlap_ids = overlap_ids[overlap_ids[:, 0] * overlap_ids[:, 1] <= 0]

    # Flip overlaps, such that the 1st column contains ids from the 1st set,
    # the 2nd column contains ids from the 2nd set.
    overlap_ids.sort(axis=-1)

    # Restore original indexes,
    overlap_ids[:, 0] = overlap_ids[:, 0] * (-1) - 1
    overlap_ids[:, 1] = overlap_ids[:, 1] - 1

    # Sort overlaps according to the 1st
    if sort:
        overlap_ids = overlap_ids[np.lexsort([overlap_ids[:, 1], overlap_ids[:, 0]])]

    return overlap_ids


def overlap_intervals(starts1, ends1, starts2, ends2, closed=False, sort=False):
    """
    Take two sets of intervals and return the indices of pairs of overlapping intervals.

    Parameters
    ----------
    starts1, ends1, starts2, ends2 : numpy.ndarray
        Interval coordinates. Warning: if provided as pandas.Series, indices
        will be ignored.

    closed : bool
        If True, then treat intervals as closed and report single-point overlaps.
    Returns
    -------
    overlap_ids : numpy.ndarray
        An Nx2 array containing the indices of pairs of overlapping intervals.
        The 1st column contains ids from the 1st set, the 2nd column has ids
        from the 2nd set.

    """

    for vec in [starts1, ends1, starts2, ends2]:
        if issubclass(type(vec), pd.core.series.Series):
            warnings.warn(
                "One of the inputs is provided as pandas.Series and its index "
                "will be ignored.",
                SyntaxWarning,
            )

    starts1 = np.asarray(starts1)
    ends1 = np.asarray(ends1)
    starts2 = np.asarray(starts2)
    ends2 = np.asarray(ends2)

    # Concatenate intervals lists
    n1 = len(starts1)
    n2 = len(starts2)
    ids1 = np.arange(0, n1)
    ids2 = np.arange(0, n2)

    # Sort all intervals together
    order1 = np.lexsort([ends1, starts1])
    order2 = np.lexsort([ends2, starts2])
    starts1, ends1, ids1 = starts1[order1], ends1[order1], ids1[order1]
    starts2, ends2, ids2 = starts2[order2], ends2[order2], ids2[order2]

    # Find interval overlaps
    match_2in1_starts = np.searchsorted(starts2, starts1, "left")
    match_2in1_ends = np.searchsorted(starts2, ends1, "right" if closed else "left")
    # "right" is intentional here to avoid duplication
    match_1in2_starts = np.searchsorted(starts1, starts2, "right")
    match_1in2_ends = np.searchsorted(starts1, ends2, "right" if closed else "left")

    # Ignore self-overlaps
    match_2in1_mask = match_2in1_ends > match_2in1_starts
    match_1in2_mask = match_1in2_ends > match_1in2_starts
    match_2in1_starts, match_2in1_ends = (
        match_2in1_starts[match_2in1_mask],
        match_2in1_ends[match_2in1_mask],
    )
    match_1in2_starts, match_1in2_ends = (
        match_1in2_starts[match_1in2_mask],
        match_1in2_ends[match_1in2_mask],
    )

    # Generate IDs of pairs of overlapping intervals
    overlap_ids = np.block(
        [
            [
                np.repeat(ids1[match_2in1_mask], match_2in1_ends - match_2in1_starts)[
                    :, None
                ],
                ids2[arange_multi(match_2in1_starts, match_2in1_ends)][:, None],
            ],
            [
                ids1[arange_multi(match_1in2_starts, match_1in2_ends)][:, None],
                np.repeat(ids2[match_1in2_mask], match_1in2_ends - match_1in2_starts)[
                    :, None
                ],
            ],
        ]
    )

    if sort:
        # Sort overlaps according to the 1st
        overlap_ids = overlap_ids[np.lexsort([overlap_ids[:, 1], overlap_ids[:, 0]])]

    return overlap_ids


def overlap_intervals_outer(starts1, ends1, starts2, ends2, closed=False):
    """
    Take two sets of intervals and return the indices of pairs of overlapping intervals,
    as well as the indices of the intervals that do not overlap any other interval.

    Parameters
    ----------
    starts1, ends1, starts2, ends2 : numpy.ndarray
        Interval coordinates. Warning: if provided as pandas.Series, indices
        will be ignored.

    closed : bool
        If True, then treat intervals as closed and report single-point overlaps.

    Returns
    -------
    overlap_ids : numpy.ndarray
        An Nx2 array containing the indices of pairs of overlapping intervals.
        The 1st column contains ids from the 1st set, the 2nd column has ids
        from the 2nd set.

    no_overlap_ids1, no_overlap_ids2 : numpy.ndarray
        Two 1D arrays containing the indices of intervals in sets 1 and 2
        respectively that do not overlap with any interval in the other set.

    """

    ovids = overlap_intervals(starts1, ends1, starts2, ends2, closed=closed)
    no_overlap_ids1 = np.where(
        np.bincount(ovids[:, 0], minlength=starts1.shape[0]) == 0
    )[0]
    no_overlap_ids2 = np.where(
        np.bincount(ovids[:, 1], minlength=starts2.shape[0]) == 0
    )[0]
    return ovids, no_overlap_ids1, no_overlap_ids2


def merge_intervals(starts, ends, min_dist=0):
    """
    Merge overlapping intervals.

    Parameters
    ----------
    starts, ends : numpy.ndarray
        Interval coordinates. Warning: if provided as pandas.Series, indices
        will be ignored.

    min_dist : float or None
        If provided, merge intervals separated by this distance or less.
        If None, do not merge non-overlapping intervals. Using
        min_dist=0 and min_dist=None will bring different results.
        bioframe uses semi-open intervals, so interval pairs [0,1) and [1,2)
        do not overlap, but are separated by a distance of 0. Such intervals
        are not merged when min_dist=None, but are merged when min_dist=0.

    Returns
    -------
    cluster_ids : numpy.ndarray
        The indices of interval clusters that each interval belongs to.
    cluster_starts : numpy.ndarray
    cluster_ends : numpy.ndarray
        The spans of the merged intervals.

    Notes
    -----
    From
    https://stackoverflow.com/questions/43600878/merging-overlapping-intervals/58976449#58976449
    """

    for vec in [starts, ends]:
        if issubclass(type(vec), pd.core.series.Series):
            warnings.warn(
                "One of the inputs is provided as pandas.Series and its index "
                "will be ignored.",
                SyntaxWarning,
            )

    starts = np.asarray(starts)
    ends = np.asarray(ends)

    order = np.lexsort([ends, starts])
    starts, ends = starts[order], ends[order]

    ends = np.maximum.accumulate(ends)
    cluster_borders = np.zeros(len(starts) + 1, dtype=bool)
    cluster_borders[0] = True
    cluster_borders[-1] = True

    if min_dist is not None:
        cluster_borders[1:-1] = starts[1:] > ends[:-1] + min_dist
    else:
        cluster_borders[1:-1] = starts[1:] >= ends[:-1]

    cluster_ids_sorted = np.cumsum(cluster_borders)[:-1] - 1
    cluster_ids = np.full(starts.shape[0], -1)
    cluster_ids[order] = cluster_ids_sorted

    cluster_starts = starts[:][cluster_borders[:-1]]
    cluster_ends = ends[:][cluster_borders[1:]]

    return cluster_ids, cluster_starts, cluster_ends


def complement_intervals(
    starts,
    ends,
    bounds=(0, np.iinfo(np.int64).max),
):

    _, merged_starts, merged_ends = merge_intervals(starts, ends, min_dist=0)

    lo = np.searchsorted(merged_ends, bounds[0], "right")
    hi = np.searchsorted(merged_starts, bounds[1], "left")

    merged_starts = merged_starts[lo:hi]
    merged_ends = merged_ends[lo:hi]

    # Trim the complement to the bounds.
    complement_starts = np.r_[bounds[0], merged_ends]
    complement_ends = np.r_[merged_starts, bounds[1]]
    lo = 1 if (complement_starts[0] >= complement_ends[0]) else 0
    hi = -1 if (complement_starts[-1] >= complement_ends[-1]) else None
    complement_starts = complement_starts[lo:hi]
    complement_ends = complement_ends[lo:hi]

    return complement_starts, complement_ends


def _closest_intervals_nooverlap(
    starts1, ends1, starts2, ends2, direction, tie_arr=None, k=1
):
    """
    For every interval in set 1, return the indices of k closest intervals
    from set 2 to the left from the interval (with smaller coordinate).
    Overlapping intervals from set 2 are not reported, unless they overlap by a single point.

    Parameters
    ----------
    starts1, ends1, starts2, ends2 : numpy.ndarray
        Interval coordinates. Warning: if provided as pandas.Series, indices
        will be ignored.

    direction : str ("left" or "right")
        Orientation of closest interval search

    tie_arr : numpy.ndarray or None
        Extra data describing intervals in set 2 to break ties when multiple intervals
        are located at the same distance. An interval with the *lowest* value is
        selected.

    k : int
        The number of neighbors to report.

    Returns
    -------
    ids: numpy.ndarray
        One Nx2 array containing the indices of pairs of closest intervals,
        reported for the neighbors in specified direction (by genomic coordinate).
        The two columns are the inteval ids from set 1, ids of the closest intevals from set 2.

    """

    for vec in [starts1, ends1, starts2, ends2]:
        if issubclass(type(vec), pd.core.series.Series):
            warnings.warn(
                "One of the inputs is provided as pandas.Series "
                "and its index will be ignored.",
                SyntaxWarning,
            )

    starts1 = np.asarray(starts1)
    ends1 = np.asarray(ends1)
    starts2 = np.asarray(starts2)
    ends2 = np.asarray(ends2)

    n1 = starts1.shape[0]
    n2 = starts2.shape[0]

    ids = np.zeros((0, 2), dtype=int)

    if k > 0 and direction=="left":
        if tie_arr is None:
            ends2_sort_order = np.argsort(ends2)
        else:
            ends2_sort_order = np.lexsort([-tie_arr, ends2])

        ids2_endsorted = np.arange(0, n2)[ends2_sort_order]
        ends2_sorted = ends2[ends2_sort_order]

        left_closest_endidx = np.searchsorted(ends2_sorted, starts1, "right")
        left_closest_startidx = np.maximum(left_closest_endidx - k, 0)

        int1_ids = np.repeat(
            np.arange(n1), left_closest_endidx - left_closest_startidx
        )
        int2_sorted_ids = arange_multi(
            left_closest_startidx, left_closest_endidx
        )

        ids = np.vstack(
            [
                int1_ids,
                ids2_endsorted[int2_sorted_ids],
                # ends2_sorted[int2_sorted_ids] - starts1[int1_ids],
                # arange_multi(left_closest_startidx - left_closest_endidx, 0)
            ]
        ).T

    elif k > 0 and direction=="right":
        if tie_arr is None:
            starts2_sort_order = np.argsort(starts2)
        else:
            starts2_sort_order = np.lexsort([tie_arr, starts2])

        ids2_startsorted = np.arange(0, n2)[starts2_sort_order]
        starts2_sorted = starts2[starts2_sort_order]

        right_closest_startidx = np.searchsorted(starts2_sorted, ends1, "left")
        right_closest_endidx = np.minimum(
            right_closest_startidx + k, n2
        )

        int1_ids = np.repeat(
            np.arange(n1), right_closest_endidx - right_closest_startidx
        )
        int2_sorted_ids = arange_multi(
            right_closest_startidx, right_closest_endidx
        )
        ids = np.vstack(
            [
                int1_ids,
                ids2_startsorted[int2_sorted_ids],
                #  starts2_sorted[int2_sorted_ids] - ends1[int1_ids],
                #  arange_multi(1, right_closest_endidx -
                #                  right_closest_startidx + 1)
            ]
        ).T


    return ids


def closest_intervals(
    starts1,
    ends1,
    starts2=None,
    ends2=None,
    k=1,
    tie_arr=None,
    ignore_overlaps=False,
    ignore_upstream=False,
    ignore_downstream=False,
    direction=None
):
    """
    For every interval in set 1, return the indices of k closest intervals from set 2.

    Parameters
    ----------
    starts1, ends1, starts2, ends2 : numpy.ndarray
        Interval coordinates. Warning: if provided as pandas.Series, indices
        will be ignored. If start2 and ends2 are None, find closest intervals
        within the same set.

    k : int
        The number of neighbors to report.

    tie_arr : numpy.ndarray or None
        Extra data describing intervals in set 2 to break ties when multiple intervals
        are located at the same distance. Intervals with *lower* tie_arr values will
        be given priority.

    ignore_overlaps : bool
        If True, ignore set 2 intervals that overlap with set 1 intervals.

    ignore_upstream, ignore_downstream : bool
        If True, ignore set 2 intervals upstream/downstream of set 1 intervals.

    direction : numpy.ndarray with dtype bool or None
        Strand vector to define the upstream/downstream orientation of the intervals.

    Returns
    -------
    closest_ids : numpy.ndarray
        An Nx2 array containing the indices of pairs of closest intervals.
        The 1st column contains ids from the 1st set, the 2nd column has ids
        from the 2nd set.

    """

    # Get overlapping intervals:
    if ignore_overlaps:
        overlap_ids = np.zeros((0, 2), dtype=int)
    elif (starts2 is None) and (ends2 is None):
        starts2, ends2 = starts1, ends1
        overlap_ids = overlap_intervals(starts1, ends1, starts2, ends2)
        overlap_ids = overlap_ids[overlap_ids[:, 0] != overlap_ids[:, 1]]
    else:
        overlap_ids = overlap_intervals(starts1, ends1, starts2, ends2)

    # Get non-overlapping intervals:
    n = len(starts1)
    all_ids = np.arange(n)

    # + directed intervals
    ids_left_upstream = _closest_intervals_nooverlap(
        starts1[direction],
        ends1[direction],
        starts2,
        ends2,
        direction="left",
        tie_arr=tie_arr,
        k=0 if ignore_upstream else k
    )
    ids_right_downstream = _closest_intervals_nooverlap(
        starts1[direction],
        ends1[direction],
        starts2,
        ends2,
        direction="right",
        tie_arr=tie_arr,
        k=0 if ignore_downstream else k
    )
    # - directed intervals
    ids_right_upstream = _closest_intervals_nooverlap(
        starts1[~direction],
        ends1[~direction],
        starts2,
        ends2,
        direction="right",
        tie_arr=tie_arr,
        k=0 if ignore_upstream else k
    )
    ids_left_downstream = _closest_intervals_nooverlap(
        starts1[~direction],
        ends1[~direction],
        starts2,
        ends2,
        direction="left",
        tie_arr=tie_arr,
        k=0 if ignore_downstream else k
    )

    # Reconstruct original indexes (b/c we split regions by direction above)
    ids_left_upstream[:, 0]   = all_ids[direction][ids_left_upstream[:, 0]]
    ids_right_downstream[:, 0] = all_ids[direction][ids_right_downstream[:, 0]]
    ids_left_downstream[:, 0] = all_ids[~direction][ids_left_downstream[:, 0]]
    ids_right_upstream[:, 0]   = all_ids[~direction][ids_right_upstream[:, 0]]

    left_ids = np.concatenate([ids_left_upstream, ids_left_downstream])
    right_ids = np.concatenate([ids_right_upstream, ids_right_downstream])

    # Increase the distance by 1 to distinguish between overlapping
    # and non-overlapping set 2 intervals.
    left_dists = starts1[left_ids[:, 0]] - ends2[left_ids[:, 1]] + 1
    right_dists = starts2[right_ids[:, 1]] - ends1[right_ids[:, 0]] + 1

    closest_ids = np.vstack([left_ids, right_ids, overlap_ids])
    closest_dists = np.concatenate(
        [left_dists, right_dists, np.zeros(overlap_ids.shape[0])]
    )

    # Sort by distance to set 1 intervals and, if present, by the tie-breaking
    # data array.
    if tie_arr is None:
        order = np.lexsort([closest_ids[:, 1], closest_dists, closest_ids[:, 0]])
    else:
        order = np.lexsort(
            [closest_ids[:, 1], tie_arr, closest_dists, closest_ids[:, 0]]
        )

    closest_ids = closest_ids[order, :2]

    # For each set 1 interval, select up to k closest neighbours.
    interval1_run_border_mask = closest_ids[:-1, 0] != closest_ids[1:, 0]
    interval1_run_borders = np.where(np.r_[True, interval1_run_border_mask, True])[0]
    interval1_run_starts = interval1_run_borders[:-1]
    interval1_run_ends = interval1_run_borders[1:]
    closest_ids = closest_ids[
        arange_multi(
            interval1_run_starts,
            lengths=np.minimum(k, interval1_run_ends - interval1_run_starts),
        )
    ]

    return closest_ids


def coverage_intervals_rle(starts, ends, weights=None):
    n = starts.shape[0]

    if weights is None:
        weights = np.ones(n, dtype=np.int64)

    borders = np.r_[starts, ends]
    coverage_change = np.r_[weights, -1 * weights]

    borders_order = np.argsort(borders)
    borders = borders[borders_order]
    coverage = np.cumsum(coverage_change[borders_order])

    return borders, coverage


def stack_intervals(starts, ends):
    n = starts.shape[0]

    borders = np.r_[starts, ends]
    lens = np.r_[ends - starts, ends - starts]
    border_types = np.r_[np.ones_like(starts), -1 * np.ones_like(ends)]
    border_ids = np.r_[np.arange(1, n + 1), -1 * np.arange(1, n + 1)]

    border_order = np.lexsort([-lens, border_types, borders])

    borders, border_ids = borders[border_order], border_ids[border_order]

    occupancy = np.zeros(2, dtype=bool)
    levels = -1 * np.ones(n, dtype=np.int64)
    for border, border_id in zip(borders, border_ids):
        interval_id = np.abs(border_id) - 1
        if border_id > 0:
            if occupancy.sum() == occupancy.shape[0]:
                occupancy = np.r_[occupancy, np.zeros_like(occupancy)]
            new_level = np.where(occupancy == False)[0][0]
            levels[interval_id] = new_level
            occupancy[new_level] = True
        if border_id < 0:
            occupancy[levels[interval_id]] = False

    return levels
