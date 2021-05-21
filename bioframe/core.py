import pandas as pd
import numpy as np
import collections
from . import ops
from .region import from_, _get_default_colnames, add_UCSC_name_column


def _verify_columns(df, colnames, return_as_bool=False):
    """
    Raises ValueError if columns with colnames are not present in dataframe df.

    df: pandas.DataFrame

    colnames: list of columns

    return_as_bool : bool
        If true, returns as a boolean instead of raising errors. Default False.

    """

    if not type(df) is pd.core.frame.DataFrame:
        if return_as_bool:
            return False
        raise ValueError("df is not a dataframe")

    if not set(colnames).issubset(df.columns):
        if return_as_bool:
            return False
        raise ValueError(
            ", ".join(set(colnames).difference(set(df.columns)))
            + " not in keys of df.columns"
        )
    if return_as_bool:
        return True


def _verify_column_dtypes(df, cols=None, return_as_bool=False):
    """
    Checks that dataframe `df` has chrom, start, end columns with valid dtypes.
    Raises TypeErrors if cols have invalid dtypes.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    return_as_bool : bool
        If true, returns as a boolean instead of raising errors. Default False.


    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    if not _verify_columns(df, [ck1, sk1, ek1], return_as_bool=True):
        if return_as_bool:
            return False
        raise ValueError("could not verify columns")

    chrom_dtype, start_dtype, end_dtype = df.dtypes[[ck1, sk1, ek1]]

    if not is_chrom_dtype(chrom_dtype):
        if return_as_bool:
            return False
        raise TypeError(
            "invalid df['chrom'] dtype, must be object, string, or categorical"
        )
    if not pd.api.types.is_integer_dtype(start_dtype):
        if return_as_bool:
            return False
        raise TypeError("invalid df['start'] dtype, must be integer")

    if not pd.api.types.is_integer_dtype(end_dtype):
        if return_as_bool:
            return False
        raise TypeError("invalid df['end'] dtype, must be integer")

    if return_as_bool:
        return True


def is_chrom_dtype(chrom_dtype):
    """
    Returns True if dtype is any of the allowed bioframe chrom dtypes, False otherwise.
    """
    return np.any(
        [
            pd.api.types.is_string_dtype(chrom_dtype),
            pd.api.types.is_object_dtype(chrom_dtype),
            pd.api.types.is_categorical_dtype(chrom_dtype),
        ]
    )


def is_bedframe(
    df,
    raise_errors=False,
    cols=None,
):
    """
    Checks that a genomic interval dataframe `df` has:
    - chrom, start, end columns
    - columns have valid dtypes (object/string/categorical, int, int)
    - all starts < ends.

    raise_errors:bool
        If true, raises errors instead of returning a boolean False for invalid properties.
        Default false.

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
            raise TypeError("Invalid column names")
        return False

    if not _verify_column_dtypes(df, cols=[ck1, sk1, ek1], return_as_bool=True):
        if raise_errors:
            raise TypeError("Invalid column dtypes")
        return False

    if ((df[ek1] - df[sk1]) < 0).any():
        if raise_errors:
            raise ValueError(
                "Invalid genomic interval dataframe: starts exceed ends for "
                + str(np.sum(((df[ek1] - df[sk1]) < 0)))
                + " intervals"
            )
        return False

    return True


def is_overlapping(df, cols=None):
    """
    tests if any genomic intervals in a bioframe `df` overlap

    Returns
    -------
    is_overlapping:bool

    """

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    df_merged = ops.merge(df, cols=cols)

    total_interval_len = np.sum((df[ek1] - df[sk1]).values)
    total_interval_len_merged = np.sum((df_merged[ek1] - df_merged[sk1]).values)

    if total_interval_len > total_interval_len_merged:
        return True
    else:
        return False


def is_viewframe(region_df, raise_errors=False, view_name_col="name", cols=None):
    """
    Checks that region_df is a valid view, namely:
    - it satisfies requirements for a bedframe, including columns for ('chrom', 'start', 'end')
    - it has an additional column, view_name_col, with default 'name'
    - it does not contain null values
    - entries in the view_name_col are unique.
    - intervals are non-overlapping

    raise_errors:bool
        If true, raises errors instead of returning a boolean for invalid properties.
        Default false.

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


def make_viewframe(
    regions,
    check_bounds=None,
    view_names_as_UCSC=False,
    view_name_col="name",
    cols=None,
):
    """
    Makes and validates a dataframe view_df, where supported input types for regions are:
    - a dictionary where keys are strings and values are integers {str:int},
    specifying regions (chrom, 0, end, chrom)
    - a pandas series of chromosomes lengths with index specifying region names
    - a list of tuples [(chrom,start,end), ...] or [(chrom,start,end,name), ...]
    - a pandas DataFrame, skips to validation step

    check_bounds : None, or chromosome sizes provided as any of valid formats above
        Optional, if provided checks if regions in the view are contained by regions
        supplied in check_bounds, typically provided as a series of chromosome sizes.
        Default None.

    view_name_col:str
        Specifies column name of the view regions. Default 'name'.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    view_df:dataframe satisfying properties of a view

    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    view_df = from_(regions, name_col=view_name_col, cols=cols)

    if check_bounds is not None:
        bounds_df = from_(check_bounds, name_col="bounds", cols=cols)
        if not is_contained(
            view_df,
            bounds_df,
            df_view_col=view_name_col,
            view_name_col="bounds",
            cols=cols,
        ):
            raise ValueError(
                "Invalid input to make a viewFrame, regions not contained by bounds"
            )

    if view_names_as_UCSC:
        view_df = add_UCSC_name_column(view_df, name_col=view_name_col, cols=cols)

    if is_viewframe(view_df, view_name_col=view_name_col, cols=cols, raise_errors=True):
        return view_df
    else:
        raise ValueError("could not make valid viewFrame, retry with new input")


def is_cataloged(
    df, view_df, raise_errors=False, df_view_col="view_region", view_name_col="name"
):
    """
    tests if all regions names in a bioframe `df` are present in the view `view_df`.

    df : pandas.DataFrame

    view_df : pandas.DataFrame

    df_view_col: str
        Name of column from df that indicates region in view.

    view_name_col: str
        Name of column from view that specifies unique region name.

    Returns
    -------
    is_cataloged:bool

    """
    if not _verify_columns(df, [df_view_col], return_as_bool=True):
        if raise_errors is True:
            raise ValueError(f"Could not find ‘{df_view_col}’ column in df")
        return False

    if not _verify_columns(view_df, [view_name_col], return_as_bool=True):
        if raise_errors is True:
            raise ValueError(f"Could not find ‘{view_name_col}’ column in view_df")
        return False

    if not set(df[df_view_col].values).issubset(set(view_df[view_name_col].values)):
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


def is_contained(
    df,
    view_df,
    raise_errors=False,
    df_view_col="view_region",
    view_name_col="name",
    cols=None,
):
    """
    tests if all genomic intervals in a bioframe `df` are cataloged and do not extend beyond their
    associated region in the view `view_df`.

    df : pandas.DataFrame

    view_df : pandas.DataFrame
        Valid viewframe.

    df_view_col:
        Column from df used to associate interviews with view regions.
        Default `view_region`.

    cols: (str, str, str)
        Column names for chrom, start, end in df.

    Returns
    -------
    is_contained:bool

    """

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    if not is_cataloged(
        df, view_df, df_view_col=df_view_col, view_name_col=view_name_col
    ):
        if raise_errors:
            raise ValueError("df not cataloged in view_df")
        return False

    df_trim = ops.trim(
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
    tests if a view `view_df` is covered by the set of genomic intervals in the bedframe `df`
    this is true if the complement is empty.

    Note this does not depend on regions assigned to intervals in df, if any, since regions are re-assigned in complement.

    Returns
    -------
    is_covering:bool

    """

    if ops.complement(
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
    tests if a view `view_df` is tiled by the set of genomic intervals in the bedframe `df`
    this is true if:
    - df is not overlapping
    - df is covering view_df
    - df is contained in view_df

    Returns
    -------
    is_tiling:bool

    """

    view_df = make_viewframe(view_df, view_name_col=view_name_col, cols=cols)

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


def sort_bedframe(
    df,
    view_df=None,
    infer_assignment=True,
    reset_index=True,
    df_view_col="view_region",
    view_name_col="name",
    cols=None,
):
    """
    Sorts a bedframe df. If no view_df is provided, sorts by cols.
    If view_df is provided, sorts df[df_view_col] by view_df[view_name_col].
    If view_df is provided but a column matching df_view_col is not in df, attempts
    to assign interavls to the view region with the largest overlap and then sorts.

    df : pandas.DataFrame
        Valid bedframe.

    view_df : pandas.DataFrame
        Valid input to make a viewframe.

    infer_assignment : bool
        If True, tries to assign df intervals to the view region with the largest overlap.
        Default True.

    reset_index : bool
        Default True.

    df_view_col:
        Column from df used to associate interviews with view regions.
        Default `view_region`.

    view_name_col:
        Column from view_df with names of regions.
        Default `name`.

    Returns
    -------
    out_df : sorted bedframe

    Notes
    -------
        df_view_col is currently returned as an ordered categorical

    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    if not is_bedframe(df, cols=cols):
        raise ValueError("not a valid bedframe, cannot sort")

    out_df = df.copy()
    if view_df is None:
        out_df.sort_values([ck1, sk1, ek1], inplace=True)

    else:
        view_df = make_viewframe(view_df, view_name_col=view_name_col, cols=cols)

        if infer_assignment and (
            not _verify_columns(out_df, [df_view_col], return_as_bool=True)
        ):
            out_df = assign_view(
                out_df,
                view_df,
                df_view_col=df_view_col,
                view_name_col=view_name_col,
                cols=cols,
            )

        if not _verify_columns(out_df, [df_view_col], return_as_bool=True):
            raise ValueError("no df_view_col not present in df, cannot sort by view")

        if not is_cataloged(
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

    if reset_index:
        out_df.reset_index(inplace=True, drop=True)

    return out_df


def is_sorted(
    df,
    view_df=None,
    infer_assignment=True,
    reset_index=True,
    df_view_col="view_region",
    view_name_col="name",
    cols=None,
):
    """
    Tests if a bedframe is changed by sorting.

    Returns
    -------
    is_sorted : bool

    """

    df_sorted = sort_bedframe(
        df.copy(),
        view_df=view_df,
        infer_assignment=infer_assignment,
        reset_index=reset_index,
        df_view_col=df_view_col,
        view_name_col=view_name_col,
        cols=cols,
    )

    if df.equals(df_sorted):
        return True
    else:
        return False


def sanitize_bedframe(
    df1,
    recast_dtypes=True,
    drop_null=False,
    start_exceed_end_action=None,
    cols=None,
):
    """
    Attempts to clean a genomic interval dataframe to be valid.

    drop_null : bool
        Drops rows with pd.NA. Default False.

    start_exceed_end_action : str or None
        Options: 'flip' or 'drop' or None. Default None.
        If 'flip', attempts to sanitize by flipping intervals with start>end.
        If 'drop' attempts to sanitize dropping intervals with start>end.
        If None, does not alter these intervals if present.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Notes
    -------
        The option start_exceed_end_action='flip' may be useful for gff files with strand information but starts > ends.

    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    out_df = df1.copy()

    _verify_columns(out_df, [ck1, sk1, ek1])

    if recast_dtypes:
        chrom_dtype, start_dtype, end_dtype = out_df.dtypes[[ck1, sk1, ek1]]
        if not is_chrom_dtype(chrom_dtype):
            out_df[ck1] = out_df[ck1].astype(str)
        if not ((start_dtype is pd.Int64Dtype()) and (end_dtype is pd.Int64Dtype())):
            out_df[sk1] = out_df[sk1].astype(pd.Int64Dtype())
            out_df[ek1] = out_df[ek1].astype(pd.Int64Dtype())

    if drop_null:
        out_df = out_df.iloc[pd.isna(out_df).any(axis=1).values == 0, :]
        out_df.reset_index(drop=True, inplace=True)

    if start_exceed_end_action is not None:
        start_exceed_end_action = start_exceed_end_action.lower()
        if ((out_df[ek1] - out_df[sk1]) < 0).any():
            inds = ((out_df[ek1] - out_df[sk1]) < 0).values
            if start_exceed_end_action == "drop":
                out_df = out_df.loc[inds == 0]
            elif start_exceed_end_action == "flip":
                out_df.loc[inds, [sk1, ek1]] = out_df.loc[inds, [ek1, sk1]].values
            else:
                raise ValueError("unknown action for intervals with start>end")
            out_df.reset_index(drop=True, inplace=True)

    if is_bedframe(out_df, cols=cols):
        return out_df
    else:
        raise ValueError("could not sanitize")


def assign_view(
    df,
    view_df,
    drop_unassigned=False,
    df_view_col="view_region",
    view_name_col="name",
    cols=None,
):
    """
    Associates genomic intervals in bedframe df with regions in viewframe view_df, based on their largest overlap.
    Note this resets the index of the returned dataframe.

    drop_unassigned : bool
        If True, drop intervals in df that are do no Default False.

    Returns
    -------
    out_df:dataframe with the associated view region for each interval in out_df[view_name_col]

    """

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    df = df.copy()
    df.reset_index(inplace=True, drop=True)

    is_bedframe(df, raise_errors=True, cols=cols)
    view_df = make_viewframe(view_df, view_name_col=view_name_col, cols=cols)

    overlap_view = ops.overlap(
        df,
        view_df,
        how="left",
        suffixes=("", "_view"),
        return_overlap=True,
        keep_order=True,
        return_index=True,
        cols1=cols,
        cols2=cols,
    )

    overlap_columns = overlap_view.columns
    overlap_view["overlap_length"] = (
        overlap_view["overlap_" + ek1] - overlap_view["overlap_" + sk1]
    )

    out_df = (
        overlap_view.sort_values("overlap_length", ascending=False)
        .drop_duplicates("index", keep="first")
        .sort_index()
    )

    out_df.rename(columns={view_name_col + "_view": df_view_col}, inplace=True)

    if drop_unassigned:
        out_df = out_df.iloc[pd.isna(out_df).any(axis=1).values == 0, :]
    out_df.reset_index(inplace=True, drop=True)

    return_cols = list(df.columns) + [df_view_col]

    return out_df[return_cols]
