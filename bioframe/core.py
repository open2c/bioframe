import pandas as pd
import numpy as np
import collections
from . import ops

_rc = {"colnames": {"chrom": "chrom", "start": "start", "end": "end"}}


def _get_default_colnames():
    return _rc["colnames"]["chrom"], _rc["colnames"]["start"], _rc["colnames"]["end"]


class update_default_colnames:
    def __init__(self, new_colnames):
        self._old_colnames = dict(_rc["colnames"])
        if isinstance(new_colnames, collections.Iterable):
            if len(new_colnames) != 3:
                raise ValueError(
                    "Please, specify new columns using a list of "
                    "3 strings or a dict!"
                )
            (
                _rc["colnames"]["chrom"],
                _rc["colnames"]["start"],
                _rc["colnames"]["end"],
            ) = new_colnames
        elif isinstance(new_colnames, collections.Mapping):
            _rc["colnames"].update(
                {
                    k: v
                    for k, v in new_colnames.items()
                    if k in ["chrom", "start", "end"]
                }
            )
        else:
            raise ValueError(
                "Please, specify new columns using a list of " "3 strings or a dict!"
            )

    def __enter__(self):
        return self

    def __exit__(self, *args):
        _rc["colnames"] = self._old_colnames


def _verify_columns(df, colnames, return_as_bool=False):
    """
    df: pandas.DataFrame

    colnames: list of columns
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
            raise TypeError("Invalid columns")
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
    Checks that region_df is a valid view, namely that it:
    - satisfies requirements for a bedframe, including columns for ('chrom', 'start', 'end')
    - has an additional column providing a unique name, view_name_col, with default 'name'
    - does not contain null values
    - that entries in the view_name_col are unique.
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
            raise TypeError("Invalid view: column names cannot be verified")
        return False

    if not is_bedframe(region_df, cols=cols):
        if raise_errors:
            raise ValueError("Invalid view: bedframe properties cannot be verified")
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
    regions, infer_chroms_from_regions=True, view_name_col="name", cols=None
):
    """
    Makes and validates a dataframe view_df, where supported input types for regions are:
    - a dictionary where keys are strings and values are integers {str:int},
    specifying regions (chrom, 0, end, name)
    - a dictionary where keys are strings and values are tuples of integers {str:(int,int)},
    specifying regions (chrom, start, end, name), e.g. {'chr1':10, 'chr2':20} or {'chr1':(0,10),'chr2':(0,20)}.
    - a pandas series of chromosomes lengths with index specifying region names
    - a dataFrame, skips to validation step

    infer_chroms_from_regions:bool
        attemps to strip 'p' or 'q' from chrom string. if False, view_name_col specifies chrom as well.
        default True.

    view_name_col:str
        specifies

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    view_df:dataframe satisfying properties of a view

    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    if type(regions) is dict:
        data = []
        for k, v in dict(regions).items():
            name = k
            if infer_chroms_from_regions:
                chrom = k.split("_")[0].replace("p", "").replace("q", "")
            else:
                chrom = k
            if isinstance(v, (tuple, list, np.ndarray)):
                start = v[0]
                end = v[1]
            elif np.isscalar(v):
                start = 0
                end = v
            else:
                raise ValueError("Unknown dict format: {type(v)}")
            data.append([chrom, start, end, name])

        view_df = pd.DataFrame(data, columns=[ck1, sk1, ek1, view_name_col])

    elif type(regions) is pd.core.series.Series:
        chroms = regions.index.values
        if infer_chroms_from_regions:
            chroms = [i.split("_")[0].replace("p", "").replace("q", "") for i in chroms]
        data = {
            ck1: chroms,
            sk1: 0,
            ek1: regions.values,
            view_name_col: regions.index.values,
        }
        view_df = pd.DataFrame(data)

    elif type(regions) is pd.core.frame.DataFrame:
        view_df = regions.copy()

    else:
        raise ValueError("Unknown region type: {type(v)}")

    if is_viewframe(
        view_df, view_name_col=view_name_col, cols=(ck1, sk1, ek1), raise_errors=True
    ):
        return view_df
    else:
        raise ValueError("could not make valid viewFrame, retry with new input")


def is_cataloged(
    df, view_df, raise_errors=False, df_view_col="view_region", view_name_col="name"
):
    """
    tests if the set of regions names in a bioframe `df` are cataloged in the view `view_df`.

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
            raise ValueError("df_view_col could not be verified in df")
        return False
    if not _verify_columns(view_df, [view_name_col], return_as_bool=True):
        if raise_errors is True:
            raise ValueError("view_name_col could not be verified in view_df")
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
    tests if all genomic intervals in a bioframe `df` are contained in the view `view_df`.

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


def is_tiling(df, view_df, df_view_col="view_region", view_name_col="name", cols=None):
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
        print("overlaps")
        return False
    if not is_covering(df, view_df, view_name_col=view_name_col, cols=None):
        print("not covered")
        return False
    if not is_contained(
        df, view_df, df_view_col=df_view_col, view_name_col=view_name_col, cols=None
    ):
        print("not contained")
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
        Valid input for a viewframe or.

    infer_assignment : bool

    reset_index : bool

    df_view_col:
        Column from df used to associate interviews with view regions.
        Default `view_region`.

    view_name_col:
        Column from view_df with names of regions.
        Default `name`.

    Returns
    -------
    out_df : sorted bedframe

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
                cols=None,
            )

        if not is_cataloged(
            out_df, view_df, df_view_col=df_view_col, view_name_col=view_name_col
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
    drop_null=True,
    flip_invalid=True,
    cols=None,
):
    """
    Attempts to clean a genomic interval dataframe to be valid.

    drop_null : bool
        Drops rows with pd.NA. Default True.

    drop_invalid:bool
        Drops invalid intervals from returned bedframe. Default False.

    flip_invalid:bool
        Flips intervals where start<end in returned bedframe.
        Default True.

    ### TODO: 
    - discuss dropping of intervals based on view_df & df_view_col
    - whether to include return_sorted, and if so, how.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

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

    if flip_invalid:
        if ((out_df[ek1] - out_df[sk1]) < 0).any():
            inds = ((out_df[ek1] - out_df[sk1]) < 0).values
            out_df.loc[inds, [sk1, ek1]] = out_df.loc[inds, [ek1, sk1]].values

    if is_bedframe(out_df, cols=cols):
        return out_df
    else:
        raise ValueError("could not sanitize")


def assign_view(
    df, view_df, df_view_col="view_region", view_name_col="name", cols=None
):
    """
    Associates genomic intervals in bedframe df with regions in viewframe view_df, based on their largest overlap.
    Currently drops all genomic intervals that are not cataloged in view_df.

    ## TODO: discuss if we should add a drop_unassigned=True option.

    Returns
    -------
    out_df:dataframe with view region in view_name_col

    """

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    is_bedframe(df, raise_errors=True, cols=cols)
    view_df = make_viewframe(view_df, view_name_col=view_name_col, cols=cols)

    overlap_idxs = ops.overlap(
        df,
        view_df,
        how="left",
        suffixes=("", "_view"),
        return_overlap=True,
        return_index=True,
        keep_order=True,
    )
    overlap_idxs["overlap"] = (
        overlap_idxs["overlap_end"] - overlap_idxs["overlap_start"]
    )

    out_df = (
        overlap_idxs.groupby("index", sort=False).apply(
            lambda group: group.nlargest(1, columns="overlap")
        )
    ).reset_index(drop=True)
    out_df.rename(columns={view_name_col + "_view": df_view_col}, inplace=True)

    return_cols = list(df.columns) + [df_view_col]
    return out_df[return_cols]
