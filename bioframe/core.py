import pandas as pd
import numpy as np
import collections
from bioframe.ops import _get_default_colnames, _verify_columns
import bioframe.definitions

def _verify_column_dtypes(df, cols=None):
    """
    Checks that a dataframe `df` has chrom, start, end columns and with valid dtypes.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    """

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    _verify_columns(df, [ck1, sk1, ek1])

    column_dtypes = df.dtypes
    if not np.any(
        [
            pd.api.types.is_string_dtype(column_dtypes[ck1]),
            pd.api.types.is_object_dtype(column_dtypes[ck1]),
            pd.api.types.is_categorical_dtype(column_dtypes[ck1]),
        ]
    ):
        raise TypeError(
            "invalid df['chrom'] dtype, must be object, string, or categorical"
        )
    if not pd.api.types.is_integer_dtype(column_dtypes[sk1]):
        raise TypeError("invalid df['start'] dtype, must be integer")
    if not pd.api.types.is_integer_dtype(column_dtypes[ek1]):
        raise TypeError("invalid df['end'] dtype, must be integer")

def is_valid(df, cols=None):
    """
    Checks that a genomic interval dataframe `df` has:
    - chrom, start, end columns
    - columns have valid dtypes
    - all starts < ends.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.
    """

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    _verify_column_dtypes(df)

    if ((df[ek1] - df[sk1]) < 0).any():
        raise ValueError("Invalid genomic interval dataframe: "+
            " starts exceed ends for "+ str(np.sum(((df[ek1] - df[sk1]) < 0)))+ " intervals")

    return True

def sanitize(df, return_sorted=True, drop_invalid=True, flip_invalid=True, cols=None):
    """
    Attempts to clean a genomic interval dataframe to be valid. 

    <<<TODO: discuss if this should accept a tuple of column names instead>>>
    return_sorted:bool
        returns the dataframe sorted by (chrom,start,end)

    drop_invalid:bool
        returns a copy of the dataframe without invalid invervals. default False.

    flip_invalid:bool
        flips intervals where start<end and returns a copy of the dataframe where these intervals have been flipped.
        default False.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    _verify_columns(df, [ck1, sk1, ek1])

    out_df = df.copy()
    if drop_invalid:
        out_df = out_df.iloc[((df[ek1] - df[sk1]) >= 0).values]

    if ((df[ek1] - df[sk1]) < 0).any() and flip_invalid:
        inds = ((df[ek1] - df[sk1]) < 0).values
        out_df.loc[inds, [sk1, ek1]] = out_df.loc[inds, [ek1, sk1]].values

    if is_valid(out_df):
        return out_df
    else:
        raise ValueError("could not sanitize")

def is_view(region_df, region_name_col="name", cols=None):
    """
    Checks that a region_df:
    - has chrom, start, end, region_name_col
    - that entries in the region_name_col are unique.
    - intervals are non-overlapping


    """

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    _verify_columns(region_df, [ck1, sk1, ek1, region_name_col])

    is_valid(region_df)

    if pd.isna(region_df).values.any():
        raise ValueError("Invalid region dataframe: cannot contain NAs")

    if len(set(region_df[region_name_col])) < len(region_df[region_name_col].values):
        raise ValueError("Invalid region dataframe: entries in region_df[region_name_col] must be unique")

    if bioframe.definitions.is_overlapping(region_df):
        raise ValueError("Invalid region dataframe: entries must be non-overlapping")

    return True

def make_view(regions, infer_chroms_from_regions=True, region_name_col="name", cols=None):
    """
    Makes and validates a dataframe view_df, where supported input types for regions are:
    - a dictionary where keys are strings and values are integers {str:int},
    specifying regions (chrom, 0, end, name)
    - a dictionary where keys are strings and values are tuples of integers {str:(int,int)},
    specifying regions (chrom, start, end, name)
    - a pandas series of chromosomes lengths with index specifying region names
    - a dataFrame, skips to validation step

    infer_chroms_from_regions:bool
        attemps to strip 'p' or 'q' from chrom string. if False, region_name_col specifies chrom as well.
        default True.

    region_name_col:str
        specifies

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.


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

        view_df = pd.DataFrame(data, columns=[ck1, sk1, ek1, region_name_col])

    elif type(regions) is pd.core.series.Series:
        chroms = regions.index.values
        if infer_chroms_from_regions:
            chroms = [i.split("_")[0].replace("p", "").replace("q", "") for i in chroms]
        data = {ck1: chroms,
                sk1: 0,
                ek1: regions.values,
                region_name_col: regions.index.values}
        view_df = pd.DataFrame(data)

    elif type(regions) is pd.core.frame.DataFrame:
        view_df = regions.copy()

    else:
        raise ValueError("Unknown region type: {type(v)}")

    is_view(view_df, region_name_col=region_name_col, cols=(ck1, sk1, ek1))

    return view_df
