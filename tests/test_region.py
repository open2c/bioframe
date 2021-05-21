from bioframe.region import parse_region
from bioframe.region import from_
import pandas as pd
import pytest


def test_parse_region():
    # UCSC-style names
    assert parse_region("chr21") == ("chr21", 0, None)
    assert parse_region("chr21:1000-2000") == ("chr21", 1000, 2000)
    assert parse_region("chr21:1,000-2,000") == ("chr21", 1000, 2000)

    # Ensembl style names
    assert parse_region("6") == ("6", 0, None)
    assert parse_region("6:1000-2000") == ("6", 1000, 2000)
    assert parse_region("6:1,000-2,000") == ("6", 1000, 2000)

    # FASTA style names
    assert parse_region("gb|accession|locus") == ("gb|accession|locus", 0, None)
    assert parse_region("gb|accession|locus:1000-2000") == (
        "gb|accession|locus",
        1000,
        2000,
    )
    assert parse_region("gb|accession|locus:1,000-2,000") == (
        "gb|accession|locus",
        1000,
        2000,
    )

    # Punctuation in names (aside from :)
    assert parse_region("name-with-hyphens-") == ("name-with-hyphens-", 0, None)
    assert parse_region("GL000207.1") == ("GL000207.1", 0, None)
    assert parse_region("GL000207.1:1000-2000") == ("GL000207.1", 1000, 2000)

    # Trailing dash
    assert parse_region("chr21:1000-") == ("chr21", 1000, None)

    # Humanized units
    assert parse_region("6:1kb-2kb") == ("6", 1000, 2000)
    assert parse_region("6:1k-2000") == ("6", 1000, 2000)
    assert parse_region("6:1kb-2M") == ("6", 1000, 2000000)
    assert parse_region("6:1Gb-") == ("6", 1000000000, None)

    with pytest.raises(ValueError):
        parse_region("chr1:2,000-1,000")  # reverse selection

    with pytest.raises(ValueError):
        parse_region("chr1::1000-2000")  # more than one colon


def test_from_():

    ### tests copied from old parse_regions
    # main functionality: convert to dataframe and create name
    df = pd.DataFrame(
        {"chrom": [f"chr{i}" for i in range(3)], "start": [1, 2, 3], "end": [4, 5, 6]}
    )
    parsed = from_(df, names_as_UCSC=True, name_col="regions")
    assert parsed.iloc[0]["chrom"] == "chr0"
    assert parsed.iloc[0]["regions"] == "chr0:1-4"

    # re-create dataframe from UCSC name alone
    df2 = pd.DataFrame({"regions": parsed["regions"].values})
    assert (from_(df2, names_as_UCSC=True, name_col="regions") == parsed).all().all()

    # re-parsing results yields the same
    assert (from_(parsed) == parsed).all().all()

    # None or False will be parsed
    assert from_([("chr1", None, 5)], fill_null={"chr1": 10})["start"].values[0] == 0

    # pull end from chromsizes
    p2 = from_([("chr1", 5, None)], fill_null={"chr1": 40}, names_as_UCSC=True)
    assert list(p2.values[0]) == ["chr1", 5, 40, "chr1:5-40"]

    # We could keep things as None if chromsizes were not proviced
    p3 = from_(["chr1", "chr2"], fill_null=False)
    assert list(p3.values[0]) == ["chr1", None, None, "chr1"]

    # we can force CUSC names
    p4 = from_([("chr1", 0, 5)], names_as_UCSC=True)
    assert list(p4.values[0]) == ["chr1", 0, 5, "chr1:0-5"]

    # nothing happens: name was autocompleted
    p5 = from_([("chr1", 0, 5)], names_as_UCSC=False)
    assert list(p5.values[0]) == ["chr1", 0, 5, "chr1"]

    # forcing UCSC names
    p7 = from_([("chr1")], fill_null={"chr1": 5}, names_as_UCSC=True)
    assert list(p7.values[0]) == ["chr1", 0, 5, "chr1:0-5"]

    # kept the strange name
    p8 = from_(["chr1:1,000,000-4M"])
    assert list(p8.values[0]) == ["chr1", 1000000, 4000000, "chr1:1,000,000-4M"]

    p9 = from_(["chr1"])
    assert list(p9.values[0]) == ["chr1", None, None, "chr1"]

    with pytest.raises(ValueError):
        from_([("ch1", 1, 2, "chr1:1-2", "puppies")])  # puppies are not allowed

    with pytest.raises(ValueError):
        from_([("chr1", 5, None)], fill_null={"chr2": 40})
