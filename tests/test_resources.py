import pandas as pd

import bioframe


def test_fetch_chromsizes():
    db = "hg38"
    for provider in ["local", "ucsc"]:
        chromsizes = bioframe.fetch_chromsizes(db, provider=provider)
        assert isinstance(chromsizes, pd.Series)
        assert chromsizes.name == "length"
        assert len(chromsizes) == 25

        chromsizes_df = bioframe.fetch_chromsizes(db, provider=provider, as_bed=True)
        assert isinstance(chromsizes_df, pd.DataFrame)
        assert list(chromsizes_df.columns) == ["chrom", "start", "end"]
        assert len(chromsizes_df) == 25

    # Check synonymous local assemblies
    assert bioframe.fetch_chromsizes("hg38", provider="local").equals(
        bioframe.fetch_chromsizes("GRCh38", provider="local")
    )


def test_fetch_chromsizes_local_vs_ucsc():
    for db in ["hg19", "hg38", "mm9", "mm10"]:
        assert bioframe.fetch_chromsizes(db, provider="local").equals(
            bioframe.fetch_chromsizes(db, provider="ucsc")
        )


def test_fetch_centromeres():
    for db in ["hg19", "hg38"]:
        # Note: UCSC will usually have a different ordering of chromosomes
        for provider in ["local", "ucsc"]:
            centromeres = bioframe.fetch_centromeres(db, provider=provider)
            assert isinstance(centromeres, pd.DataFrame)
            assert list(centromeres.columns) == ["chrom", "start", "end", "mid"]
            assert len(centromeres) == 24
