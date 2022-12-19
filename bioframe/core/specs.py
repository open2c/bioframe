import collections
import pandas as pd
import numpy as np

__all__ = [
    "update_default_colnames",
    "is_chrom_dtype",
]

_rc = {"colnames": {"chrom": "chrom", "start": "start", "end": "end"}}


def _get_default_colnames():
    """
    Returns default column names.

    These defaults be updated with :func:`update_default_colnames`.

    Returns
    -------
    colnames : triplet (str, str, str)

    """
    return _rc["colnames"]["chrom"], _rc["colnames"]["start"], _rc["colnames"]["end"]


class update_default_colnames:
    def __init__(self, new_colnames):
        self._old_colnames = dict(_rc["colnames"])
        if isinstance(new_colnames, collections.abc.Iterable):
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
        elif isinstance(new_colnames, collections.abc.Mapping):
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


def _verify_columns(df, colnames, unique_cols=False, return_as_bool=False):
    """
    Raises ValueError if columns with colnames are not present in dataframe df.

    Parameters
    ----------
    df: pandas.DataFrame

    colnames: list of column names

    return_as_bool : bool
        If True, returns as a boolean instead of raising errors. Default False.

    """

    if not type(df) is pd.core.frame.DataFrame:
        if return_as_bool:
            return False
        raise ValueError("df is not a dataframe")

    if unique_cols:
        if len(set(colnames)) < len(colnames):
            raise ValueError("column names must be unique")

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

    Parameters
    ----------
    df : pandas.DataFrame

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
