from functools import partial
from urllib.parse import urljoin, urlencode
import urllib
import os
import posixpath as pp
import os.path as op
import numpy as np
import pandas as pd
import requests
import socket
import base64
import glob

import pkg_resources

from .schemas import SCHEMAS, UCSC_MRNA_FIELDS
from .fileops import (
    read_table,
    read_chromsizes,
)

__all__ = [
    "fetch_chromsizes",
    "fetch_centromeres",
    "UCSCClient",
]


def _check_connectivity(reference="http://www.google.com"):
    try:
        urllib.request.urlopen(reference, timeout=5)
        return True
    except urllib.request.URLError:
        return False
    except socket.timeout:
        return False


def fetch_chromsizes(
    db,
    provider=None,
    filter_chroms=True,
    chrom_patterns=(r"^chr[0-9]+$", r"^chr[XY]$", r"^chrM$"),
    natsort=True,
    as_bed=False,
    **kwargs,
):
    """
    Fetch chromsizes from the UCSC database or local storage.

    Parameters
    ----------
    provider : str
        The provider of chromsizes. Currently, only 'ucsc' is implemented.
    filter_chroms : bool, optional
        Filter for chromosome names given in ``chrom_patterns``.
    chrom_patterns : sequence, optional
        Sequence of regular expressions to capture desired sequence names.
    natsort : bool, optional
        Sort each captured group of names in natural order. Default is True.
    as_bed : bool, optional
        If True, return chromsizes as an interval dataframe (chrom, start, end).
    **kwargs :
        Passed to :func:`pandas.read_csv`

    Returns
    -------
    Series of integer bp lengths indexed by sequence name or an interval dataframe.

    """

    if provider == "local":
        fpath = f"data/{db}.chrom.sizes"
        if pkg_resources.resource_exists("bioframe.io", fpath):
            return read_chromsizes(
                pkg_resources.resource_filename("bioframe.io", fpath)
            )
        else:
            raise LookupError(f"Assembly '{db}' not found in local storage")

    if provider == "ucsc" or provider is None:
        return UCSCClient(db).fetch_chromsizes(
            filter_chroms=filter_chroms,
            chrom_patterns=chrom_patterns,
            natsort=natsort,
            as_bed=as_bed,
            **kwargs,
        )
    else:
        raise ValueError("Unknown provider '{}'".format(provider))


def _origins_from_cytoband(cyb, band_col="gieStain"):
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
            raise ValueError(
                f"Expected 2 'acen' bands for {chrom}, found {len(group)}"
            )
        acens = group.sort_values("start")
        cens.append({
            "chrom": chrom,
            "start": acens.iloc[0]["start"],
            "end": acens.iloc[1]["end"],
            "mid": acens.iloc[0]["end"],
        })
    return pd.DataFrame.from_records(cens)


def _origins_from_ucsccentromeres(cens):
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
    cens = cens.groupby("chrom").agg({
        "start": np.min, 
        "end": np.max
    }).reset_index()
    cens["mid"] = (cens["start"] + cens["end"]) // 2
    cens = (
        cens[["chrom", "start", "end", "mid"]]
        .sort_values("chrom")
        .reset_index(drop=True)
    )
    return cens


def fetch_centromeres(db, provider=None, merge=True, verbose=False):
    """
    Extract centromere locations for a given assembly 'db' from a variety
    of file formats in UCSC (cytoband, centromeres) depending on
    availability, returning a DataFrame.

    Parameters
    ----------

    db : str

    merge : bool
        Whether to merge all centromere intervals per chromosome into
        one consolidated centromere interval.
        Default True.

    Returns
    -------
    DataFrame with centromere 'chrom', 'start', 'end', 'mid'.

    Notes
    -----
    The priority goes as
    - Local (not implemented)
    - centromeres.txt
    - cytoBandIdeo
    - cytoBand
    
    Gap files no longer provide centromere information.
    Currently only works for human assemblies.

    """
    if provider == "local":
        raise NotImplementedError("local method not currently implemented")
        fpath = f"data/{db}.centromeres"
        if pkg_resources.resource_exists("bioframe.io", fpath):
            return read_chromsizes(
                pkg_resources.resource_filename("bioframe.io", fpath)
            )
        else:
            raise LookupError(f"Centromeres for '{db}' not found in local storage")

    if provider == "ucsc" or provider is None:
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
            raise ValueError("No source for centromere data found.")
    else:
        raise NotImplementedError("currently UCSC is only implemented provider")
    
    if schema == "centromeres":
        return _origins_from_ucsccentromeres(df)
    else:
        return _origins_from_cytoband(df)


class UCSCClient:
    BASE_URL = "http://hgdownload.cse.ucsc.edu/"

    def __init__(self, db):
        self._db = db
        self._db_url = urljoin(self.BASE_URL, f"goldenPath/{db}/")

    def fetch_chromsizes(
        self,
        filter_chroms=True,
        chrom_patterns=(r"^chr[0-9]+$", r"^chr[XY]$", r"^chrM$"),
        natsort=True,
        as_bed=False,
        **kwargs,
    ):
        url = urljoin(self._db_url, f"bigZips/{self._db}.chrom.sizes")
        return read_chromsizes(
            url,
            filter_chroms=filter_chroms,
            chrom_patterns=chrom_patterns,
            natsort=natsort,
            as_bed=as_bed,
            **kwargs,
        )

    def fetch_centromeres(self, **kwargs):
        url = urljoin(self._db_url, "database/centromeres.txt.gz")
        return read_table(url, schema="centromeres")

    def fetch_gaps(self, **kwargs):
        url = urljoin(self._db_url, "database/gap.txt.gz")
        return read_table(
            url,
            schema="gap",
            usecols=["chrom", "start", "end", "length", "type", "bridge"],
            **kwargs,
        )

    def fetch_cytoband(self, ideo=False, **kwargs):
        if ideo:
            url = urljoin(self._db_url, "database/cytoBandIdeo.txt.gz")
        else:
            url = urljoin(self._db_url, "database/cytoBand.txt.gz")
        return read_table(url, schema="cytoband")
    
    def fetch_mrna(self, **kwargs):
        url = urljoin(self._db_url, "database/all_mrna.txt.gz")
        return read_table(
            url,
            schema=UCSC_MRNA_FIELDS,
            **kwargs,
        )


