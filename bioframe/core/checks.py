import pandas as pd
import numpy as np
from . import construction
from .specs import _get_default_colnames, _verify_columns, _verify_column_dtypes
from .. import ops

__all__ = [
    "is_bedframe",
    "is_cataloged",
    "is_overlapping",
    "is_viewframe",
    "is_contained",
    "is_covering",
    "is_tiling",
    "is_sorted",
]


def is_bedframe(
    df,
    raise_errors=False,
    cols=None,
):
    """
    Checks that required bedframe properties are satisfied for dataframe `df`.

    This includes:

    - chrom, start, end columns
    - columns have valid dtypes (object/string/categorical, int/pd.Int64Dtype, int/pd.Int64Dtype)
    - for each interval, if any of chrom, start, end are null, then all are null
    - all starts < ends.

    Parameters
    ----------
    df : pandas.DataFrame

    raise_errors : bool
        If True, raises errors instead of returning a boolean False for invalid properties.
        Default False.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    is_bedframe:bool

    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    if not _verify_columns(df, [ck1, sk1, ek1], return_as_bool=True):
        if raise_errors:
            raise TypeError("Invalid bedFrame: Invalid column names")
        return False

    if not _verify_column_dtypes(df, cols=[ck1, sk1, ek1], return_as_bool=True):
        if raise_errors:
            raise TypeError("Invalid bedFrame: Invalid column dtypes")
        return False

    nan_intervals = pd.isnull(df[[ck1, sk1, ek1]])
    if (~(~nan_intervals.any(axis=1) | nan_intervals.all(axis=1))).any():
        if raise_errors:
            raise ValueError(
                "Invalid bedFrame: Invalid null values (if any of chrom, start, end are null, then each must be null)"
            )
        return False

    if ((df[ek1] - df[sk1]) < 0).any():
        if raise_errors:
            raise ValueError(
                "Invalid bedFrame: starts exceed ends for "
                + str(np.sum(((df[ek1] - df[sk1]) < 0)))
                + " intervals"
            )
        return False

    return True


def is_cataloged(
    df, view_df, raise_errors=False, df_view_col="view_region", view_name_col="name"
):
    """
    Tests if all region names in `df[df_view_col]` are present in `view_df[view_name_col]`.

    Parameters
    ----------
    df : pandas.DataFrame

    view_df : pandas.DataFrame

    raise_errors : bool
        If True, raises errors instead of returning a boolean False for invalid properties.
        Default False.

    df_view_col: str
        Name of column from df that indicates region in view.

    view_name_col: str
        Name of column from view that specifies  region name.

    Returns
    -------
    is_cataloged:bool

    Notes
    -----
    Does not check if names in `view_df[view_name_col]` are unique.

    """
    if not _verify_columns(df, [df_view_col], return_as_bool=True):
        if raise_errors is True:
            raise ValueError(f"Could not find ‘{df_view_col}’ column in df")
        return False

    if not _verify_columns(view_df, [view_name_col], return_as_bool=True):
        if raise_errors is True:
            raise ValueError(f"Could not find ‘{view_name_col}’ column in view_df")
        return False

    if not set(df[df_view_col].copy().dropna().values).issubset(
        set(view_df[view_name_col].values)
    ):
        if raise_errors is True:
            raise ValueError(
                "The following regions in df[df_view_col] not in view_df[view_name_col]: \n"
                + "{}".format(
                    set(df[df_view_col].values).difference(
                        set(view_df[view_name_col].values)
                    )
                )
            )
        return False

    return True


def is_overlapping(df, cols=None):
    """
    Tests if any genomic intervals in a bioframe `df` overlap.

    Also see :func:`bioframe.ops.merge()`.

    Parameters
    ----------
    df : pandas.DataFrame

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    is_overlapping:bool

    """
    from ..ops import merge

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    df_merged = merge(df, cols=cols)

    total_interval_len = np.sum((df[ek1] - df[sk1]).values)
    total_interval_len_merged = np.sum((df_merged[ek1] - df_merged[sk1]).values)

    if total_interval_len > total_interval_len_merged:
        return True
    else:
        return False


def is_viewframe(region_df, raise_errors=False, view_name_col="name", cols=None):
    """
    Checks that `region_df` is a valid viewFrame.

    This includes:

    - it satisfies requirements for a bedframe, including columns for ('chrom', 'start', 'end')
    - it has an additional column, view_name_col, with default 'name'
    - it does not contain null values
    - entries in the view_name_col are unique.
    - intervals are non-overlapping

    Parameters
    ----------

    region_df : pandas.DataFrame
        Dataframe of genomic intervals to be tested.

    raise_errors : bool
        If True, raises errors instead of returning a boolean False for invalid properties.
        Default False.

    view_name_col : str
        Specifies column name of the view regions. Default 'name'.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    is_viewframe:bool

    """

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    if not _verify_columns(
        region_df, [ck1, sk1, ek1, view_name_col], return_as_bool=True
    ):
        if raise_errors:
            raise TypeError("Invalid view: invalid column names")
        return False

    if not is_bedframe(region_df, cols=cols):
        if raise_errors:
            raise ValueError("Invalid view: not a bedframe")
        return False

    if pd.isna(region_df).values.any():
        if raise_errors:
            raise ValueError("Invalid view: cannot contain NAs")
        return False

    if len(set(region_df[view_name_col])) < len(region_df[view_name_col].values):
        if raise_errors:
            raise ValueError(
                "Invalid view: entries in region_df[view_name_col] must be unique"
            )
        return False

    if is_overlapping(region_df, cols=cols):
        if raise_errors:
            raise ValueError("Invalid view: entries must be non-overlapping")
        return False

    return True


def is_contained(
    df,
    view_df,
    raise_errors=False,
    df_view_col=None,
    view_name_col="name",
    cols=None,
):
    """
    Tests if all genomic intervals in a bioframe `df` are cataloged and do not extend beyond their
    associated region in the view `view_df`.

    Parameters
    ----------
    df : pandas.DataFrame

    view_df : pandas.DataFrame
        Valid viewframe.

    raise_errors : bool
        If True, raises errors instead of returning a boolean False for invalid properties.
        Default False.

    df_view_col:
        Column from df used to associate interviews with view regions.
        Default `view_region`.

    cols: (str, str, str)
        Column names for chrom, start, end in df.

    Returns
    -------
    is_contained:bool

    """
    from ..ops import trim

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    if df_view_col is None:
        try:
            df_view_assigned = ops.overlap(df, view_df)
            assert (df_view_assigned["end_"].isna()).sum() == 0
            assert (df_view_assigned["start_"].isna()).sum() == 0
            assert (df_view_assigned["end"] <= df_view_assigned["end_"]).all()
            assert (df_view_assigned["start"] >= df_view_assigned["start_"]).all()
        except AssertionError:
            if raise_errors:
                raise AssertionError("df not contained in view_df")
            else:
                return False
        return True

    if not is_cataloged(
        df, view_df, df_view_col=df_view_col, view_name_col=view_name_col
    ):
        if raise_errors:
            raise ValueError("df not cataloged in view_df")
        return False

    df_trim = trim(
        df, view_df=view_df, df_view_col=df_view_col, view_name_col=view_name_col
    )
    is_start_trimmed = np.any(df[sk1].values != df_trim[sk1].values)
    is_end_trimmed = np.any(df[ek1].values != df_trim[ek1].values)

    if is_start_trimmed or is_end_trimmed:
        if raise_errors:
            raise ValueError("df not contained in view_df")
        return False
    else:
        return True


def is_covering(df, view_df, view_name_col="name", cols=None):
    """
    Tests if a view `view_df` is covered by the set of genomic intervals in the bedframe `df`.

    This test is true if ``complement(df,view_df)`` is empty. Also note this test ignores regions assigned to
    intervals in `df` since regions are re-assigned in :func:`bioframe.ops.complement`.

    Parameters
    ----------
    df : pandas.DataFrame

    view_df : pandas.DataFrame
        Valid viewFrame.

    view_name_col:
        Column from view_df with view region names. Default `name`.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    is_covering:bool

    """
    from ..ops import complement

    if complement(
        df,
        view_df=view_df,
        view_name_col=view_name_col,
        cols=cols,
    ).empty:
        return True
    else:
        return False


def is_tiling(
    df,
    view_df,
    raise_errors=False,
    df_view_col="view_region",
    view_name_col="name",
    cols=None,
):
    """
    Tests if a view `view_df` is tiled by the set of genomic intervals in the bedframe `df`.

    This is true if:

    - df is not overlapping
    - df is covering view_df
    - df is contained in view_df

    Parameters
    ----------
    df : pandas.DataFrame

    view_df : pandas.DataFrame
        valid viewFrame

    raise_errors : bool
        If True, raises errors instead of returning a boolean False for invalid properties.
        Default False.

    df_view_col: str
        Name of column from df that indicates region in view.

    view_name_col: str
        Name of column from view that specifies unique region name.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    is_tiling:bool

    """

    view_df = construction.make_viewframe(
        view_df, view_name_col=view_name_col, cols=cols
    )

    if is_overlapping(df):
        if raise_errors:
            raise ValueError("overlaps")
        return False
    if not is_covering(df, view_df, view_name_col=view_name_col, cols=None):
        if raise_errors:
            raise ValueError("not covered")
        return False
    if not is_contained(
        df, view_df, df_view_col=df_view_col, view_name_col=view_name_col, cols=None
    ):
        if raise_errors:
            raise ValueError("not contained")
        return False
    return True


def is_sorted(
    df,
    view_df=None,
    reset_index=True,
    df_view_col=None,
    view_name_col="name",
    cols=None,
):
    """
    Tests if a bedframe is changed by sorting.

    Also see :func:`bioframe.ops.sort_bedframe`.

    Parameters
    ----------
    df : pandas.DataFrame

    view_df : pandas.DataFrame | dict-like
        Optional view to pass to ``sort_bedframe``.
        When it is dict-like :func:'bioframe.make_viewframe' will
        be used to convert to viewframe. If view_df is not provided
        df is assumed to be sorted by chrom and start.

    reset_index : bool
        Optional argument to pass to ``sort_bedframe``.

    df_view_col: None | str
        Name of column from df that indicates region in view.
        If None, :func:'bioframe.assign_view' will be used to assign view regions.
        Default None.

    view_name_col: str
        Name of column from view that specifies unique region name.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    is_sorted : bool

    """
    from ..ops import sort_bedframe

    df_sorted = sort_bedframe(
        df.copy(),
        view_df=view_df,
        reset_index=reset_index,
        df_view_col=df_view_col,
        view_name_col=view_name_col,
        cols=cols,
    )

    if df.equals(df_sorted):
        return True
    else:
        return False
