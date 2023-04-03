from typing import List, Tuple, Union, Literal, Dict
from dataclasses import dataclass
import pkg_resources

import pandas as pd
import numpy as np
import yaml

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


def assemblies_available() -> pd.DataFrame:
    path = pkg_resources.resource_filename("bioframe.io", "data/_assemblies.yml")
    with open(path) as f:
        assemblies = yaml.safe_load(f)
    return pd.DataFrame.from_records(assemblies)


def assembly_info(
        name: str,
        seq_types: Union[List, Tuple, Literal["all"]] = None, 
        seq_units: Union[List, Tuple, Literal["all"]] = None
    ) -> GenomeAssembly:
    """
    Get information about a genome assembly.

    Parameters
    ----------
    name : str
        Name of the assembly. If the name contains a dot, it is interpreted as
        a provider name and a version, e.g. "hg38". Otherwise, the provider
        is inferred from the default provider for the organism.
    seq_types : list or tuple or "all", optional
        Sequence types to include in the assembly. If not specified, the
        default sequence types for the assembly are used.
    seq_units : list or tuple or "all", optional
        Assembly units to include in the assembly. If not specified, only 
        sequences from the default units for the assembly are included.

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
    
    >>> assembly_info("hg38", seq_types=("assembled", "non-nuclear"))

    >>> assembly_info("ucsc.hg38", seq_units=("unplaced",))

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
    default_seq_types = assembly["default_types"]
    default_seq_units = assembly["default_units"]
    seqinfo_path = assembly["seqinfo"]
    seqinfo = pd.read_table(
        pkg_resources.resource_filename("bioframe.io", f"data/{seqinfo_path}")
    )
    mask = np.ones(len(seqinfo), dtype=bool)
    if seq_types is None:
        mask &= seqinfo["type"].isin(default_seq_types)
    elif isinstance(seq_types, (tuple, list)):
        mask &= seqinfo["type"].isin(seq_types)
    elif isinstance(seq_types, str) and seq_types != "all":
        raise ValueError(f"seq_types must be a tuple or 'all', not {seq_types}")
    if seq_units is None:
        mask &= seqinfo["unit"].isin(default_seq_units)
    elif isinstance(seq_units, (tuple, list)):
        mask &= seqinfo["unit"].isin(seq_units)
    elif isinstance(seq_units, str) and seq_units != "all":
        raise ValueError(f"seq_units must be a tuple or 'all', not {seq_units}")
    seqinfo = seqinfo.loc[mask]

    return GenomeAssembly(
        organism=assembly["organism"],
        provider=assembly["provider"],
        provider_version=assembly["provider_version"],
        release_year=assembly["release_year"],
        seqinfo=seqinfo,
        url=assembly["url"],
    )
