import pandas as pd
import pytest

from bioframe.io.assembly import GenomeAssembly, assemblies_available, assembly_info


def test_assemblies_available():
    assemblies = assemblies_available()
    assert isinstance(assemblies, pd.DataFrame)
    for col in ["provider", "provider_build", "default_roles", "default_units"]:
        assert col in assemblies.columns


def test_assembly_info():
    hg38 = assembly_info("hg38")
    assert isinstance(hg38, GenomeAssembly)
    assert hg38.provider == "ucsc"
    assert hg38.provider_build == "hg38"
    assert isinstance(hg38.chromsizes, pd.Series)
    assert isinstance(hg38.chromnames, list)
    assert isinstance(hg38.alias_dict, dict)

    assert isinstance(hg38.seqinfo, pd.DataFrame)
    for col in ["name", "length", "aliases", "role", "unit"]:
        assert col in hg38.seqinfo.columns

    assert isinstance(hg38.viewframe, pd.DataFrame)
    for col in ["chrom", "start", "end", "name"]:
        assert col in hg38.viewframe.columns

    hg38 = assembly_info("ucsc.hg38", roles=("assembled", "unlocalized"))
    assert isinstance(hg38, GenomeAssembly)

    with pytest.raises(ValueError):
        assembly_info("ncbi.hg38")  # provider-name mismatch

    assert isinstance(hg38.cytobands, pd.DataFrame)
    for col in ["chrom", "start", "end", "band", "stain"]:
        assert col in hg38.cytobands.columns

    sacCer3 = assembly_info("sacCer3")
    assert sacCer3.cytobands is None
