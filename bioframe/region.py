import re
import pandas as pd

__all__ = ["parse_region", "parse_regions"]


def atoi(s):
    return int(s.replace(",", ""))


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


def regions_add_name_column(reg_df, cols=("chrom", "start", "end", "name")):
    """
    Checks that input dataframe has "name" column
    If not, auto-creates a UCSC name

    cols are names of columns (unlikely to change)
    """
    if cols[3] in reg_df:
        return reg_df[list(cols)]
    df = reg_df.copy()
    data = zip(df[cols[0]], df[cols[1]], df[cols[2]])
    df[cols[3]] = ["{0}:{1}-{2}".format(*i) for i in data]
    return df[list(cols)]


def parse_regions(
    regions,
    chromsizes=None,
    fill_missing=True,
    check_bounds=True,
    overwrite_names=False,
    cols=("chrom", "start", "end", "name"),
):
    """
    Coerce a sequence of genomic regions to a 4-column regions dataframe.

    Input may be a DataFrame or an iterable of region strings or sequences.
    DataFrame or list/tuple input may contain a 4th "name" field. If not, the
    region names will be generated from the coordinates using UCSC-style
    notation.

    Parameters
    ---------
    regions : dataframe or iterable
        Object to convert to a region dataframe.
    chromsizes : dict-like, optional
        Chromsizes used to fill unknown end coordinates, if needed.
    fill_missing : bool, optional
        Try to replace unknown coordinates with 0 or the chromosome length.
        If True, will raise an error if chromsizes are needed but not provided.
        If False, None coordinates are not checked.
    check_bounds : bool, optional, default: True
        If True, check that start <= end. Ignored if ``fill_missing`` is False.
    overwrite_names : bool, optional
        If True, will overwrite an exiting name column with UCSC-style region
        strings.
    cols : tuple, optional
        Names of the columns (unlikely to change)

    Returns
    -------
    DataFrame with columns - chrom, start, end, name.

    Notes
    -----
    See this gist for examples
    https://gist.github.com/mimakaev/9d2eb07dc746c6010304d795c99125ed

    """
    if isinstance(regions, pd.DataFrame):
        # DataFrame input
        in_cols = regions.columns
        if all([i in in_cols for i in cols[:3]]):
            new_regions = regions.copy()

        elif "regions" in in_cols:
            parsed = [parse_region_string(i) for i in regions["regions"].values]
            new_regions = pd.DataFrame(parsed, columns=cols[:3])
            new_regions[cols[3]] = regions["regions"]
        else:
            raise ValueError(
                f"regions not found in input dataframe with columns {in_cols}"
            )
    else:
        try:
            regions = list(regions)  # it has to be converted to list
        except:
            raise ValueError("Input should be a dataframe, or iterable")

        if isinstance(regions[0], str):
            # Sequence of region strings
            parsed = [parse_region_string(i) for i in regions]
            new_regions = pd.DataFrame(parsed, columns=cols[:3])
            new_regions[cols[3]] = regions

        else:
            # Sequence of region tuples/lists
            regions = [tuple(i) for i in regions]
            ulen = list(set([len(i) for i in regions]))
            if len(ulen) != 1:
                raise ValueError(f"Different lengths of values in input data: {ulen}")
            if ulen[0] in [3, 4]:
                new_regions = pd.DataFrame(regions, columns=cols[: ulen[0]])
            else:
                raise ValueError(f"Wrong number of columns:{ulen[0]}")

    # chrom
    new_regions[cols[0]] = new_regions[cols[0]].astype(str)

    # start, end
    if fill_missing:
        starts = []
        for i in new_regions[cols[1]].values:
            try:
                starts.append(int(i))
            except:
                if not i:
                    starts.append(0)
                else:
                    raise ValueError(f"Wrong start {i}; False or None are accepted")

        new_regions[cols[1]] = starts
        ends = []
        ends_orig = new_regions[cols[2]].values
        chroms = new_regions[cols[0]].values
        for i in range(len(new_regions)):
            try:
                ends.append(int(ends_orig[i]))
            except (TypeError, ValueError):
                if ends_orig[i]:
                    raise ValueError(
                        f"Wrong end {ends_orig[i]}; False or None are accepted"
                    )
                if chromsizes is not None:
                    if chroms[i] in chromsizes:
                        ends.append(chromsizes[chroms[i]])
                    else:
                        raise ValueError(f"{chroms[i]} not found in chromsizes")
                else:
                    raise ValueError(
                        f"End {ends_orig[i]} undefined and no chromsizes proviced"
                    )
        new_regions[cols[2]] = ends

    if check_bounds and fill_missing:
        if (new_regions[cols[2]] < new_regions[cols[1]]).any():
            raise ValueError("start > end detected")

    # name
    if cols[3] not in new_regions or overwrite_names:
        try:
            new_regions.pop("name")
        except KeyError:
            pass
        new_regions = regions_add_name_column(new_regions)
    else:
        new_regions = new_regions[list(cols)].copy()

    return new_regions
