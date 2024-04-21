import re
from typing import Optional, Tuple, Union

import pandas as pd

__all__ = [
    "parse_region",
    "parse_region_string",
    "is_complete_ucsc_string",
    "to_ucsc_string",
]

NUMERIC_REGEX = re.compile("([0-9,.]+)")

RANGE_TOKEN_SPEC = [
    ("HYPHEN", r"-"),
    ("COORD", r"[0-9,]+(\.[0-9]*)?(?:[a-z]+)?"),
    ("OTHER", r".+"),
]

RANGE_REGEX = re.compile(
    r"\s*" + r"|\s*".join(rf"(?P<{name}>{token})" for name, token in RANGE_TOKEN_SPEC),
    re.IGNORECASE,
)


def to_ucsc_string(grange: Tuple[str, int, int]) -> str:
    """
    Convert a grange to a UCSC string.

    Parameters
    ----------
    grange : tuple or other iterable
        chrom, start, end

    Returns
    -------
    str
        UCSC-style genomic range string, '{chrom}:{start}-{end}'
    """
    return "{}:{}-{}".format(*grange)


def is_complete_ucsc_string(s: str) -> bool:
    """
    Returns True if a string can be parsed into a completely informative
    (chrom, start, end) format.

    Parameters
    ----------
    s : str

    Returns
    -------
    bool
        True if able to be parsed and ``end`` is known.

    """
    if not isinstance(s, str):
        return False
    _, _, end = parse_region_string(s)
    if end is None:
        return False
    return True


def _parse_humanized_int(s: str) -> int:
    _, value, unit = NUMERIC_REGEX.split(s.replace(",", ""))

    # No multiplier unit, just return the integer value
    if not len(unit):
        return int(value)

    # Parse and apply the multiplier. Remaining decimal places are dropped.
    value = float(value)
    unit = unit.upper().strip()
    if unit in ("K", "KB"):
        value *= 1_000
    elif unit in ("M", "MB"):
        value *= 1_000_000
    elif unit in ("G", "GB"):
        value *= 1_000_000_000
    else:
        raise ValueError(f"Unknown unit '{unit}'")
    return int(value)


def parse_region_string(s: str) -> Tuple[str, int, int]:
    """
    Parse a UCSC-style genomic range string into a triple.

    Parameters
    ----------
    s : str
        UCSC-style genomic range string, e.g. "chr5:10,100,000-30,000,000".

    Returns
    -------
    tuple
        (str, int or None, int or None)

    See also
    --------
    parse_region
    """

    def _tokenize(s):
        for match in RANGE_REGEX.finditer(s):
            name = match.lastgroup
            yield name, match.group(name)

    def _parse_range(token_stream):
        name, token = next(token_stream, (None, None))
        if name != "COORD":
            raise ValueError(f"Expected COORD; got unexpected token: {name}: {token}")
        start = _parse_humanized_int(token)

        name, token = next(token_stream, (None, None))
        if name != "HYPHEN":
            raise ValueError(f"Expected HYPHEN; got unexpected token: {name}: {token}")

        name, token = next(token_stream, (None, None))
        if name is None:  # No end coordinate
            end = None
        elif name == "COORD":
            end = _parse_humanized_int(token)
        else:
            raise ValueError(f"Expected COORD; got unexpected token: {name}: {token}")

        return start, end

    parts = s.split(":")

    chrom = parts[0].strip()
    if not len(chrom):
        raise ValueError("Chromosome name cannot be empty")

    if len(parts) < 2:
        return (chrom, None, None)

    start, end = _parse_range(_tokenize(parts[1]))

    return chrom, start, end


def _parse_region_record(grange: tuple) -> Tuple[str, int, int]:
    """
    Coerce a genomic range record into a triple.

    Parameters
    ----------
    grange : str or tuple
        * A triple (chrom, start, end), where ``start`` or ``end`` may be
          ``None``.
        * A quadruple or higher-order tuple, e.g. (chrom, start, end, name).
          ``name`` and other fields will be ignored.

    Returns
    -------
    tuple
        A well-formed genomic range triple (str, int, int).
    """
    if len(grange) < 3:
        raise ValueError("Length of a range record should be at least 3")
    chrom, start, end = grange[:3]
    chrom = str(chrom)
    start = int(start) if start is not None else start
    end = int(end) if end is not None else end
    return chrom, start, end


def parse_region(
    grange: Union[str, tuple],
    chromsizes: Optional[Union[dict, pd.Series]] = None,
    *,
    check_bounds: bool = True,
) -> Tuple[str, int, int]:
    """
    Coerce a genomic range string or sequence type into a triple.

    Parameters
    ----------
    grange : str or tuple
        * A UCSC-style genomic range string, e.g. "chr5:10,100,000-30,000,000".
        * A triple (chrom, start, end), where ``start`` or ``end`` may be
          ``None``.
        * A quadruple or higher-order tuple, e.g. (chrom, start, end, name).
          ``name`` and other fields will be ignored.

    chromsizes : dict or Series, optional
        Lookup table of sequence lengths for bounds checking and for
        filling in a missing end coordinate.

    check_bounds : bool, optional [default: True]
        If True, check that the genomic range is within the bounds of the
        sequence.

    Returns
    -------
    tuple
        A well-formed genomic range triple (str, int, int).

    Notes
    -----
    Genomic ranges are interpreted as half-open intervals (0-based starts,
    1-based ends) along the length coordinate of a sequence.

    Sequence names may contain any character except for whitespace and colon.

    The start coordinate should be 0 or greater and the end coordinate should
    be less than or equal to the length of the sequence, if the latter is
    known. These are enforced when ``check_bounds`` is ``True``.

    If the start coordinate is missing, it is assumed to be 0. If the end
    coordinate is missing and chromsizes are provided, it is replaced with the
    length of the sequence.

    The end coordinate **must** be greater than or equal to the start.

    The start and end coordinates may be suffixed with k(b), M(b), or G(b)
    multipliers, case-insentive. e.g. "chr1:1K-2M" is equivalent to
    "chr1:1000-2000000".
    """
    if isinstance(grange, str):
        chrom, start, end = parse_region_string(grange)
    else:
        chrom, start, end = _parse_region_record(grange)

    # Fill in missing end coordinate if possible
    clen = None
    if chromsizes is not None:
        try:
            clen = chromsizes[chrom]
        except KeyError:
            raise ValueError(f"Unknown sequence label: {chrom}")
        if end is None:
            end = clen

    # Fill in missing start coordinate
    if start is None:
        start = 0

    if end is not None and (end < start):
        raise ValueError("End cannot be less than start")

    if check_bounds and (start < 0 or (clen is not None and end > clen)):
        raise ValueError(f"Genomic range out of bounds: [{start}, {end})")

    return chrom, start, end
