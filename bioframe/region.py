import re
import pandas as pd
import numpy as np

__all__ = ["parse_region"]


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


def atoi(s):
    return int(s.replace(",", ""))


### manipulating UCSC strings ###

def parse_humanized(s):
    _NUMERIC_RE = re.compile("([0-9,.]+)")
    _, value, unit = _NUMERIC_RE.split(s.replace(",", ""))
    if not len(unit):
        return int(value)

    value = float(value)
    unit = unit.upper().strip()
    if unit in ("K", "KB"):
        value *= 1000
    elif unit in ("M", "MB"):
        value *= 1000000
    elif unit in ("G", "GB"):
        value *= 1000000000
    else:
        raise ValueError("Unknown unit '{}'".format(unit))
    return int(value)


def parse_region_string(s):
    """
    Parse a UCSC-style genomic region string into a triple.

    Parameters
    ----------
    s : str
        UCSC-style string, e.g. "chr5:10,100,000-30,000,000". Ensembl and FASTA
        style sequence names are allowed. End coordinate must be greater than or
        equal to start.

    Returns
    -------
    (str, int or None, int or None)

    """

    def _tokenize(s):
        token_spec = [
            ("HYPHEN", r"-"),
            ("COORD", r"[0-9,]+(\.[0-9]*)?(?:[a-z]+)?"),
            ("OTHER", r".+"),
        ]
        tok_regex = r"\s*" + r"|\s*".join(r"(?P<%s>%s)" % pair for pair in token_spec)
        tok_regex = re.compile(tok_regex, re.IGNORECASE)
        for match in tok_regex.finditer(s):
            typ = match.lastgroup
            yield typ, match.group(typ)

    def _check_token(typ, token, expected):
        if typ is None:
            raise ValueError("Expected {} token missing".format(" or ".join(expected)))
        else:
            if typ not in expected:
                raise ValueError('Unexpected token "{}"'.format(token))

    def _expect(tokens):
        typ, token = next(tokens, (None, None))
        _check_token(typ, token, ["COORD"])
        start = parse_humanized(token)

        typ, token = next(tokens, (None, None))
        _check_token(typ, token, ["HYPHEN"])

        typ, token = next(tokens, (None, None))
        if typ is None:
            return start, None

        _check_token(typ, token, ["COORD"])
        end = parse_humanized(token)
        if end < start:
            raise ValueError("End coordinate less than start")

        return start, end

    parts = s.split(":")
    chrom = parts[0].strip()
    if not len(chrom):
        raise ValueError("Chromosome name cannot be empty")
    if len(parts) < 2:
        return (chrom, None, None)
    start, end = _expect(_tokenize(parts[1]))
    return (chrom, start, end)

def is_complete_UCSC_string(mystring):
    """
    Check if a string can be parsed into chrom,start,end format.
    """
    if type(mystring) is not str:
        return False
    if len(parse_region_string(mystring)) != 3:
        return False
    if parse_region_string(mystring)[2] is None:
        return False
    return True


def parse_region(reg, chromsizes=None):
    """
    Coerce a genomic region string or sequence into a triple (chrom, start, end).

    Genomic regions are represented as half-open intervals (0-based starts,
    1-based ends) along the length coordinate of a contig/scaffold/chromosome.

    Parameters
    ----------
    reg : str or tuple
        UCSC-style genomic region string, or
        Triple (chrom, start, end), where ``start`` or ``end`` may be ``None``.
        Quadriple (chrom, start, end, name) (name is ignored).
    chromsizes : mapping, optional
        Lookup table of scaffold lengths to check against ``chrom`` and the
        ``end`` coordinate. Required if ``end`` is not supplied.

    Returns
    -------
    A well-formed genomic region triple (str, int, int)

    """
    if isinstance(reg, str):
        chrom, start, end = parse_region_string(reg)
    else:
        if len(reg) not in [3, 4]:
            raise ValueError("length of a region should be 3 or 4")
        chrom, start, end = reg[:3]
        start = int(start) if start is not None else start
        end = int(end) if end is not None else end

    try:
        clen = chromsizes[chrom] if chromsizes is not None else None
    except KeyError:
        raise ValueError("Unknown sequence label: {}".format(chrom))

    start = 0 if start is None else start
    if end is None:
        end = clen  # if clen is None, end is None too!

    if (end is not None) and (end < start):
        raise ValueError("End cannot be less than start")

    if start < 0 or (clen is not None and end > clen):
        raise ValueError("Genomic region out of bounds: [{}, {})".format(start, end))

    return chrom, start, end


### conversions from various input formats into dataframes ###

def from_dict(regions, name_col="name", cols=None):
    """
    Makes a dataframe from a dictionary of {str,int} pairs, interpreted as chromosome names.
    Note that {str,(int,int)} dictionaries of tuples are no longer supported!
    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    data = []
    for k, v in dict(regions).items():
        chrom = k
        name = k
        if np.isscalar(v):
            start = 0
            end = v
        else:
            raise ValueError("Unsupported dict format: {type(v)}")
        data.append([chrom, start, end, name])
    return pd.DataFrame(data, columns=[ck1, sk1, ek1, name_col])


def from_series(regions, name_col="name", cols=None):
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    chroms = regions.index.values
    data = {
        ck1: chroms,
        sk1: 0,
        ek1: regions.values,
        name_col: regions.index.values,
    }
    return pd.DataFrame(data)


def from_list(regions, name_col="name", cols=None):
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    df = pd.DataFrame(regions)
    if df.shape[1] == 3:
        df.columns = [ck1, sk1, ek1]
        df[name_col] = df[ck1].values.copy()
    elif df.shape[1] == 4:
        df.columns = [ck1, sk1, ek1, name_col]
    else:
        raise ValueError('wrong number of columns for list input format')
    return df

def from_UCSC_string_list(region_list, name_col="name", cols=None):
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    parsed = [parse_region_string(i) for i in region_list]
    df = pd.DataFrame(parsed, columns=[ck1,sk1,ek1])
    df[name_col] = region_list
    return df

def from_(regions, names_as_UCSC=False, fill_null = False, name_col="name", cols=None):
    """
    Attempts to make a dataframe with columns [chr,start,end,name_col] from a variety of input types.
    Currently supported inputs:
    - dataframe
    - series of UCSC strings
    - dictionary of {str:int} key value pairs
    - pandas series where the index is interpreted as chromosomes and values are interpreted as end
    - list of tuples or lists, either [(chrom,start,end)] or [(chrom,start,end,name)]
    
    names_as_UCSC : bool
        replaces values in name_col with UCSC strings made from (chrom,start,end).
        Default False.

    fill_null : False or dictionary
        Accepts a dictionary of {str:int} pairs, interpreted as chromosome sizes. 
        Kept or backwards compatibility. Default False.

    name_col : str
        Column name. Default "name".

    cols : (str,str,str)
        Names for dataframe columns. 
        Default None sets them with get_default_colnames().
    
    Returns
    -------
    out_df:dataframe 

    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    if type(regions) is pd.core.frame.DataFrame:
        if set([ck1,sk1,ek1]).issubset(regions.columns):
            out_df = regions.copy()
        elif (len(regions[name_col].values.shape) == 1) and is_complete_UCSC_string(regions[name_col].values[0]): 
            out_df = from_UCSC_string_list(regions[name_col].values, name_col=name_col, cols=[ck1, sk1, ek1])
        else:
            raise ValueError("Unknown dataFrame format: check column names")
    
    elif type(regions) is dict:
        out_df = from_dict(regions, name_col=name_col, cols=[ck1, sk1, ek1])
    
    elif type(regions) is pd.core.series.Series:
        out_df = from_series(regions, name_col=name_col, cols=[ck1, sk1, ek1])
    
    elif type(regions) is list:
        if len(np.shape(regions))==1 and type(regions[0]) is str:
            out_df = from_UCSC_string_list(regions, name_col=name_col, cols=[ck1, sk1, ek1])
        else:
            out_df = from_list(regions, name_col=name_col, cols=[ck1, sk1, ek1])
    else:
        raise ValueError("Unknown input format: {type(regions)}")

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

    if names_as_UCSC:
        out_df = add_UCSC_name_column(out_df, name_col=name_col, cols=cols)

    return out_df

def to_UCSC_string(triplet):
    return "{0}:{1}-{2}".format(*triplet)

def add_UCSC_name_column(reg_df, name_col="name", cols=None):
    """
    Auto-creates a UCSC name chrom:start-end for each region (chrom,start,end) in reg_df, replacing name_col if it exists.
    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    df = reg_df.copy()
    data = zip(df[ck1], df[sk1], df[ek1])
    df[name_col] = [ to_UCSC_string(i) for i in data]
    return df

