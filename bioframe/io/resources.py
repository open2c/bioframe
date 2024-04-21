import urllib
from functools import partial
from typing import Union
from urllib.parse import urljoin

import numpy as np
import pandas as pd

from .assembly import assembly_info
from .fileops import read_chromsizes, read_table
from .schemas import SCHEMAS

__all__ = [
    "fetch_chromsizes",
    "fetch_centromeres",
    "UCSCClient",
]


def fetch_chromsizes(
    db: str,
    *,
    provider: str = "local",
    as_bed: bool = False,
    filter_chroms: bool = True,
    chrom_patterns: tuple = (r"^chr[0-9]+$", r"^chr[XY]$", r"^chrM$"),
    natsort: bool = True,
    **kwargs,
) -> Union[pd.Series, pd.DataFrame]:
    """
    Fetch chromsizes from local storage or the UCSC database.

    Parameters
    ----------
    db : str
        Assembly name.
    provider : str, optional [default: "local"]
        The provider of chromsizes. Either "local" for local storage or "ucsc".
    as_bed : bool, optional
        If True, return chromsizes as an interval DataFrame (chrom, start, end)
        instead of a Series.

    The remaining options only apply to provider="ucsc".

    filter_chroms : bool, optional
        Filter for chromosome names given in ``chrom_patterns``.
    chrom_patterns : sequence, optional
        Sequence of regular expressions to capture desired sequence names.
    natsort : bool, optional
        Sort each captured group of names in natural order. Default is True.
    **kwargs :
        Passed to :func:`pandas.read_csv`

    Returns
    -------
    Series of integer bp lengths indexed by sequence name or BED3 DataFrame.

    Notes
    -----
    For more fine-grained control over the chromsizes from local storage,
    use :func:`bioframe.assembly_info`.

    Examples
    --------
    >>> fetch_chromsizes("hg38")
    name
    chr1     248956422
    chr2     242193529
    chr3     198295559
    ...      ...
    chrX     156040895
    chrY      57227415
    chrM         16569
    Name: length, dtype: int64

    >>> fetch_chromsizes("hg38", as_bed=True)
            chrom      start        end
    0        chr1          0  248956422
    1        chr2          0  242193529
    2        chr3          0  198295559
    ...      ...
    21       chrX          0  156040895
    22       chrY          0   57227415
    23       chrM          0      16569

    See also
    --------
    bioframe.assembly_info
    bioframe.UCSCClient
    """
    if provider == "local":
        assembly = assembly_info(db)
        if as_bed:
            return assembly.viewframe[["chrom", "start", "end"]].copy()
        else:
            return assembly.chromsizes
    elif provider == "ucsc":
        return UCSCClient(db).fetch_chromsizes(
            filter_chroms=filter_chroms,
            chrom_patterns=chrom_patterns,
            natsort=natsort,
            as_bed=as_bed,
            **kwargs,
        )
    else:
        raise ValueError(f"Unknown provider '{provider}'")


def _origins_from_cytoband(
    cyb: pd.DataFrame, band_col: str = "gieStain"
) -> pd.DataFrame:
    """
    Extract chromosomal origin positions separating chromosome arms from
    cytological band data. Takes the cytological origin, i.e. the boundary
    between the two bands labeled 'acen'.

    Parameters
    ----------
    cyb : pandas.DataFrame
        DataFrame with cytoband data.

    Returns
    -------
    pandas.DataFrame
        A dataframe with columns 'chrom', 'start', 'end', 'mid'.
    """
    cyb = cyb[cyb[band_col] == "acen"]
    grouped = cyb.groupby("chrom", sort=False)
    cens = []
    for chrom, group in grouped:
        if not len(group) == 2:
            raise ValueError(f"Expected 2 'acen' bands for {chrom}, found {len(group)}")
        acens = group.sort_values("start")
        cens.append(
            {
                "chrom": chrom,
                "start": acens.iloc[0]["start"],
                "end": acens.iloc[1]["end"],
                "mid": acens.iloc[0]["end"],
            }
        )
    return pd.DataFrame.from_records(cens)


def _origins_from_ucsccentromeres(cens: pd.DataFrame) -> pd.DataFrame:
    """
    Extract chromosomal origin positions from UCSC centromeres.txt table
    describing centromere model sequences. Takes the midpoint of all
    modeled centromere sequences.

    Parameters
    ----------
    cens : pandas.DataFrame
        DataFrame with centromeres.txt data.

    Returns
    -------
    pandas.DataFrame
        A dataframe with columns 'chrom', 'start', 'end', 'mid'.
    """
    cens = cens.groupby("chrom").agg({"start": np.min, "end": np.max}).reset_index()
    cens["mid"] = (cens["start"] + cens["end"]) // 2
    cens = (
        cens[["chrom", "start", "end", "mid"]]
        .sort_values("chrom")
        .reset_index(drop=True)
    )
    return cens


def fetch_centromeres(db: str, provider: str = "local") -> pd.DataFrame:
    """
    Extract centromere locations for a given assembly 'db' from a variety
    of file formats in UCSC (cytoband, centromeres) depending on
    availability, returning a DataFrame.

    Parameters
    ----------
    db : str
        Assembly name.
    provider : str, optional [default: "local"]
        The provider of centromere data. Either "local" for local storage
        or "ucsc".

    Returns
    -------
    DataFrame with centromere 'chrom', 'start', 'end', 'mid'.

    Notes
    -----
    When provider="local", centromeres are derived from cytoband tables
    in local storage.

    Whe provider="ucsc", the fallback priority goes as follows:
    - UCSC cytoBand
    - UCSC cytoBandIdeo
    - UCSC centromeres.txt

    Note that UCSC "gap" files no longer provide centromere information.

    Currently only works for human assemblies.

    See also
    --------
    bioframe.assembly_info
    bioframe.UCSCClient
    """
    if provider == "local":
        assembly = assembly_info(db)
        cyb = assembly.cytobands
        if cyb is None:
            raise ValueError(
                f"No source for centromere data found from provider '{provider}'."
            )
        return _origins_from_cytoband(cyb, band_col="stain")

    elif provider == "ucsc":
        client = UCSCClient(db)
        fetchers = [
            ("cytoband", client.fetch_cytoband),
            ("cytoband", partial(client.fetch_cytoband, ideo=True)),
            ("centromeres", client.fetch_centromeres),
        ]

        for schema, fetcher in fetchers:
            try:
                df = fetcher()
                break
            except urllib.error.HTTPError:
                pass
        else:
            raise ValueError(
                f"No source for centromere data found from provider '{provider}'."
            )

        if schema == "centromeres":
            return _origins_from_ucsccentromeres(df)
        else:
            return _origins_from_cytoband(df)

    else:
        raise ValueError(f"Unknown provider '{provider}'")


class UCSCClient:
    BASE_URL = "https://hgdownload.soe.ucsc.edu/"

    def __init__(self, db: str):
        self._db = db
        self._db_url = urljoin(self.BASE_URL, f"goldenPath/{db}/")

    def fetch_chromsizes(
        self,
        filter_chroms: bool = True,
        chrom_patterns: tuple = (r"^chr[0-9]+$", r"^chr[XY]$", r"^chrM$"),
        natsort: bool = True,
        as_bed: bool = False,
        **kwargs,
    ) -> Union[pd.Series, pd.DataFrame]:
        url = urljoin(self._db_url, f"bigZips/{self._db}.chrom.sizes")
        return read_chromsizes(
            url,
            filter_chroms=filter_chroms,
            chrom_patterns=chrom_patterns,
            natsort=natsort,
            as_bed=as_bed,
            **kwargs,
        )

    def fetch_centromeres(self, **kwargs) -> pd.DataFrame:
        url = urljoin(self._db_url, "database/centromeres.txt.gz")
        return read_table(url, schema="centromeres", **kwargs)

    def fetch_gaps(self, **kwargs):
        url = urljoin(self._db_url, "database/gap.txt.gz")
        return read_table(
            url,
            schema="gap",
            usecols=["chrom", "start", "end", "length", "type", "bridge"],
            **kwargs,
        )

    def fetch_cytoband(self, ideo: bool = False, **kwargs) -> pd.DataFrame:
        if ideo:
            url = urljoin(self._db_url, "database/cytoBandIdeo.txt.gz")
        else:
            url = urljoin(self._db_url, "database/cytoBand.txt.gz")
        return read_table(url, schema="cytoband")

    def fetch_mrna(self, **kwargs) -> pd.DataFrame:
        url = urljoin(self._db_url, "database/all_mrna.txt.gz")
        return read_table(
            url,
            schema=SCHEMAS["all_mrna"],
            **kwargs,
        )
