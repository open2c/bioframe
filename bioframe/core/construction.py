import pandas as pd
import numpy as np
from .specs import _get_default_colnames, _verify_columns, is_chrom_dtype
from .stringops import parse_region_string, to_ucsc_string, is_complete_ucsc_string
from . import checks

__all__ = [
    "from_dict",
    "from_series",
    "from_list",
    "from_any",
    "make_viewframe",
    "sanitize_bedframe",
]

### conversions from various input formats into dataframes ###


def from_dict(regions, cols=None):
    """
    Makes a dataframe from a dictionary of {str,int} pairs, interpreted as chromosome names.

    Note that {str,(int,int)} dictionaries of tuples are no longer supported!

    Parameters
    ----------

    regions : dict

    name_col : str
        Default 'name'.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    df : pandas.DataFrame
    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    data = []
    for k, v in dict(regions).items():
        chrom = k
        if np.isscalar(v):
            start = 0
            end = v
        else:
            raise ValueError("Unsupported dict format: {type(v)}")
        data.append([chrom, start, end])
    return pd.DataFrame(data, columns=[ck1, sk1, ek1])


def from_series(regions, cols=None):
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    chroms = regions.index.values
    data = {ck1: chroms, sk1: 0, ek1: regions.values}
    return pd.DataFrame(data)


def from_list(regions, name_col="name", cols=None):
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    df = pd.DataFrame(regions)
    if df.shape[1] == 3:
        df.columns = [ck1, sk1, ek1]
    elif df.shape[1] == 4:
        df.columns = [ck1, sk1, ek1, name_col]
    else:
        raise ValueError("wrong number of columns for list input format")
    return df


def from_ucsc_string_list(region_list, cols=None):
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    parsed = [parse_region_string(i) for i in region_list]
    df = pd.DataFrame(parsed, columns=[ck1, sk1, ek1])
    return df


def from_any(regions, fill_null=False, name_col="name", cols=None):
    """
    Attempts to make a genomic interval dataframe with columns [chr, start, end, name_col] from a variety of input types.

    Parameters
    ----------
    regions : supported input
        Currently supported inputs:

            - dataframe
            - series of UCSC strings
            - dictionary of {str:int} key value pairs
            - pandas series where the index is interpreted as chromosomes and values are interpreted as end
            - list of tuples or lists, either [(chrom,start,end)] or [(chrom,start,end,name)]
            - tuple of tuples or lists, either [(chrom,start,end)] or [(chrom,start,end,name)]

    fill_null : False or dictionary
        Accepts a dictionary of {str:int} pairs, interpreted as chromosome sizes.
        Kept or backwards compatibility. Default False.

    name_col : str
        Column name. Only used if 4 column list is provided. Default "name".

    cols : (str,str,str)
        Names for dataframe columns.
        Default None sets them with get_default_colnames().

    Returns
    -------
    out_df:dataframe

    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    if type(regions) is pd.core.frame.DataFrame:
        if set([ck1, sk1, ek1]).issubset(regions.columns):
            out_df = regions.copy()
        elif (len(regions[name_col].values.shape) == 1) and is_complete_ucsc_string(
            regions[name_col].values[0]
        ):
            out_df = from_ucsc_string_list(
                regions[name_col].values, cols=[ck1, sk1, ek1]
            )
        else:
            raise ValueError("Unknown dataFrame format: check column names")

    elif type(regions) is dict:
        out_df = from_dict(regions, cols=[ck1, sk1, ek1])

    elif type(regions) is pd.core.series.Series:
        out_df = from_series(regions, cols=[ck1, sk1, ek1])

    elif type(regions) is tuple:
        if np.shape(regions) == (3,):
            out_df = from_list([regions], name_col=name_col, cols=[ck1, sk1, ek1])

        elif len(np.shape(regions)) == 1 and type(regions[0]) is str:
            out_df = from_ucsc_string_list(regions, cols=[ck1, sk1, ek1])
        else:
            out_df = from_list(list(regions), name_col=name_col, cols=[ck1, sk1, ek1])

    elif type(regions) is list:
        if np.shape(regions) == (3,):
            out_df = from_list([regions], name_col=name_col, cols=[ck1, sk1, ek1])
        elif len(np.shape(regions)) == 1 and type(regions[0]) is str:
            out_df = from_ucsc_string_list(regions, cols=[ck1, sk1, ek1])
        else:
            out_df = from_list(regions, name_col=name_col, cols=[ck1, sk1, ek1])
    else:
        raise ValueError(f"Unknown input format: {type(regions)}")

    if fill_null:
        try:
            out_df[sk1].fillna(0, inplace=True)
            ends = []
            for i in range(len(out_df)):
                if out_df[ek1].values[i] is None:
                    ends.append(fill_null[out_df[ck1].values[i]])
                else:
                    ends.append(out_df[ek1].values[i])
            out_df[ek1] = ends
        except:
            raise ValueError("could not fill ends with provided chromsizes")

    return out_df


def add_ucsc_name_column(reg_df, name_col="name", cols=None):
    """
    Auto-creates a UCSC name 'chrom:start-end' for each region (chrom,start,end) in reg_df.

    Replaces name_col if it exists.



    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    df = reg_df.copy()
    _verify_columns(df, [ck1, sk1, ek1])
    data = zip(df[ck1], df[sk1], df[ek1])
    df[name_col] = [to_ucsc_string(i) for i in data]
    return df


def make_viewframe(
    regions,
    check_bounds=None,
    name_style=None,
    view_name_col="name",
    cols=None,
):
    """
    Makes and validates a dataframe `view_df` out of regions.

    Parameters
    ----------
    regions : supported input type
        Currently supported input types:

            - a dictionary where keys are strings and values are integers {str:int}, specifying regions (chrom, 0, end, chrom)
            - a pandas series of chromosomes lengths with index specifying region names
            - a list of tuples [(chrom,start,end), ...] or [(chrom,start,end,name), ...]
            - a pandas DataFrame, skips to validation step

    name_style : None or "ucsc"
        If None and no column view_name_col, propagate values from cols[0]
        If "ucsc" and no column view_name_col, create UCSC style names

    check_bounds : None, or chromosome sizes provided as any of valid formats above
        Optional, if provided checks if regions in the view are contained by regions
        supplied in check_bounds, typically provided as a series of chromosome sizes.
        Default None.

    view_name_col : str
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

    view_df = from_any(regions, name_col=view_name_col, cols=cols)

    if check_bounds is not None:
        bounds_df = from_any(check_bounds, name_col="bounds", cols=cols)
        if not checks.is_contained(
            view_df,
            bounds_df,
            df_view_col=None,
            view_name_col="bounds",
            cols=cols,
        ):
            raise ValueError(
                "Invalid input to make a viewFrame, regions not contained by bounds"
            )

    if not view_name_col in view_df.columns:
        if name_style is None:
            view_df[view_name_col] = view_df[ck1].values
        elif name_style.lower() == "ucsc":
            view_df = add_ucsc_name_column(view_df, name_col=view_name_col, cols=cols)
        else:
            raise ValueError("unknown value for name_style")

    if checks.is_viewframe(
        view_df, view_name_col=view_name_col, cols=cols, raise_errors=True
    ):
        return view_df
    else:
        raise ValueError("could not make valid viewFrame, retry with new input")


def sanitize_bedframe(
    df1,
    recast_dtypes=True,
    drop_null=False,
    start_exceed_end_action=None,
    cols=None,
):
    """
    Attempts to clean a genomic interval dataframe to be a valid bedframe.

    Parameters
    ----------
    df1 : pandas.DataFrame

    recast_dtypes : bool
        Whether to attempt to recast column dtypes to pandas nullable dtypes.

    drop_null : bool
        Drops rows with pd.NA. Default False.

    start_exceed_end_action : str or None
        Options: 'flip' or 'drop' or None. Default None.

            - If 'flip', attempts to sanitize by flipping intervals with start>end.
            - If 'drop' attempts to sanitize dropping intervals with start>end.
            - If None, does not alter these intervals if present.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    Returns
    -------
    out_df : pandas.DataFrame
        Sanitized dataframe satisfying the properties of a bedframe.

    Notes
    ------
    The option ``start_exceed_end_action='flip'`` may be useful for gff files with strand information but starts > ends.

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

    nan_intervals = pd.isnull(out_df[[ck1, sk1, ek1]]).any(axis=1)
    out_df.loc[nan_intervals, [ck1, sk1, ek1]] = pd.NA
    if drop_null:
        out_df.dropna(axis=0, inplace=True)
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

    if checks.is_bedframe(out_df, cols=cols):
        return out_df
    else:
        raise ValueError("could not sanitize")
