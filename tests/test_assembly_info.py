from bioframe.io.assembly import list_assemblies, assembly_info, GenomeAssembly

import pandas as pd
import pytest


def test_list_assemblies():
    assemblies = list_assemblies()
    assert isinstance(assemblies, pd.DataFrame)
    for col in ["provider", "provider_version", "default_types", "default_units"]:
        assert col in assemblies.columns


def test_assembly_info():
    hg38 = assembly_info("hg38")
    assert isinstance(hg38, GenomeAssembly)
    assert isinstance(hg38.chromsizes, pd.Series)
    assert isinstance(hg38.seqinfo, pd.DataFrame)

    hg38 = assembly_info("ucsc.hg38", seq_types=("assembled", "non-nuclear"))
    assert isinstance(hg38, GenomeAssembly)

    with pytest.raises(ValueError):
        assembly_info("ncbi.hg38")
