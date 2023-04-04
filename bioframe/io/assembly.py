from typing import List, Tuple, Union, Dict
from dataclasses import dataclass
import pkg_resources

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

import pandas as pd
import numpy as np
import yaml
from bioframe import make_viewframe

__all__ = ["assemblies_available", "assembly_info"]


@dataclass
class GenomeAssembly:
    """
    A dataclass containing information about sequences in a genome assembly.    
    """
    organism: str
    provider: str
    provider_version: str
    release_year: str
    seqinfo: pd.DataFrame
    url: str
    alias_dict: Dict[str, str] = None

    def __post_init__(self):
        self.alias_dict = {}
        alias_lists = self.seqinfo["aliases"].str.split(',')
        names = self.seqinfo["name"]
        for aliases, name in zip(alias_lists, names):
            for alias in aliases:
                self.alias_dict[alias] = name

    @property
    def chromsizes(self) -> pd.Series:
        return self.seqinfo.set_index("name")["length"]

    @property
    def chromnames(self) -> List[str]:
        return self.seqinfo["name"].tolist()
    
    @property
    def viewframe(self) -> pd.DataFrame:
        return make_viewframe(self.chromsizes.to_dict())


def assemblies_available() -> pd.DataFrame:
    path = pkg_resources.resource_filename("bioframe.io", "data/_assemblies.yml")
    with open(path) as f:
        assemblies = yaml.safe_load(f)
    return pd.DataFrame.from_records(assemblies)


def assembly_info(
        name: str,
        roles: Union[List, Tuple, Literal["all"]] = None, 
        units: Union[List, Tuple, Literal["all"]] = None
    ) -> GenomeAssembly:
    """
    Get information about a genome assembly.

    Parameters
    ----------
    name : str
        Name of the assembly. If the name contains a dot, it is interpreted as
        a provider name and a version, e.g. "hg38". Otherwise, the provider
        is inferred if the version name is unique.
    roles : list or tuple or "all", optional
        Sequence roles to include in the assembly info. If not specified, only
        sequences with the default sequence roles for the assembly are shown.
    units : list or tuple or "all", optional
        Assembly units to include in the assembly info. If not specified, only 
        sequences from the default units for the assembly are shown.

    Returns
    -------
    GenomeAssembly
        A dataclass containing information about the assembly.

    Raises
    ------
    ValueError
        If the assembly name is not found or is not unique.
    
    Examples
    --------
    >>> hg38 = assembly_info("hg38")
    >>> hg38.chromsizes
    name
    chr1    248956422
    chr2    242193529
    chr3    198295559
    ...     ...
    
    >>> assembly_info("hg38", roles=("assembled", "non-nuclear"))

    >>> assembly_info("ucsc.hg38", units=("unplaced",))

    """
    assemblies = assemblies_available()
    if "." in name:
        provider, name = name.split(".", 1)
        provider = provider.lower()
    else:
        provider = None

    if provider is None:
        q = f"provider_version == '{name}'"
    else:
        q = f"provider == '{provider}' and provider_version == '{name}'"
    
    result = assemblies.query(q)
    if len(result) == 0:
        raise ValueError(f"Assembly not found: {name}")
    elif len(result) > 1:
        raise ValueError(f"Assembly identifer not unique: {result}")

    assembly = result.iloc[0]
    default_roles = assembly["default_roles"]
    default_units = assembly["default_units"]
    seqinfo_path = assembly["seqinfo"]
    seqinfo = pd.read_table(
        pkg_resources.resource_filename("bioframe.io", f"data/{seqinfo_path}")
    )
    mask = np.ones(len(seqinfo), dtype=bool)
    if roles is None:
        mask &= seqinfo["role"].isin(default_roles)
    elif isinstance(roles, (tuple, list)):
        mask &= seqinfo["role"].isin(roles)
    elif isinstance(roles, str) and roles != "all":
        raise ValueError(f"roles must be a tuple or 'all', not {roles}")
    if units is None:
        mask &= seqinfo["unit"].isin(default_units)
    elif isinstance(units, (tuple, list)):
        mask &= seqinfo["unit"].isin(units)
    elif isinstance(units, str) and units != "all":
        raise ValueError(f"units must be a tuple or 'all', not {units}")
    seqinfo = seqinfo.loc[mask]

    return GenomeAssembly(
        organism=assembly["organism"],
        provider=assembly["provider"],
        provider_version=assembly["provider_version"],
        release_year=assembly["release_year"],
        seqinfo=seqinfo,
        url=assembly["url"],
    )
