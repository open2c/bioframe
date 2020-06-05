from bioframe.region import parse_region, parse_regions
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
        parse_region("chr1::1000-2000")  # more than one colon


def test_parse_regions():
    # main functionality: convert to dataframe and create name
    df = pd.DataFrame(
        {"chrom": [f"chr{i}" for i in range(3)], "start": [1, 2, 3], "end": [4, 5, 6]}
    )
    parsed = parse_regions(df)
    assert parsed.iloc[0]["chrom"] == "chr0"
    assert parsed.iloc[0]["name"] == "chr0:1-4"

    # re-create dataframe from UCSC name alone
    df2 = pd.DataFrame({"regions": parsed["name"].values})
    assert (parse_regions(df2) == parsed).all().all()

    # re-parsing results yields the same
    assert (parse_regions(parsed) == parsed).all().all()

    # None or False will be parsed
    assert parse_regions([("chr1", None, 5)])["start"].values[0] == 0

    # pull end from chromsizes
    p2 = parse_regions([("chr1", 5, None)], chromsizes={"chr1": 40})
    assert list(p2.values[0]) == ["chr1", 5, 40, "chr1:5-40"]

    # We could keep things as None if chromsizes were not proviced
    p3 = parse_regions(["chr1", "chr2"], replace_None=False)
    assert list(p3.values[0]) == ["chr1", None, None, "chr1"]

    # we can force CUSC names
    p4 = parse_regions([("chr1", 0, 5)], force_UCSC_names=True)
    assert list(p4.values[0]) == ["chr1", 0, 5, "chr1:0-5"]

    # nothing happens: name was autocompleted
    p5 = parse_regions([("chr1", 0, 5)], force_UCSC_names=False)
    assert list(p5.values[0]) == ["chr1", 0, 5, "chr1:0-5"]

    # "chr1" parsed but interpreted as name
    p6 = parse_regions([("chr1")], chromsizes={"chr1": 5}, force_UCSC_names=False)
    assert list(p6.values[0]) == ["chr1", 0, 5, "chr1"]  # "chr1" interpreted as name

    # forcing UCSC names
    p7 = parse_regions([("chr1")], chromsizes={"chr1": 5}, force_UCSC_names=True)
    assert list(p7.values[0]) == ["chr1", 0, 5, "chr1:0-5"]

    # kept the strange name
    p8 = parse_regions(["chr1:1,000,000-4M"])
    assert list(p8.values[0]) == ["chr1", 1000000, 4000000, "chr1:1,000,000-4M"]

    p9 = parse_regions(["chr1"], replace_None=False)
    assert list(p9.values[0]) == ["chr1", None, None, "chr1"]

    with pytest.raises(ValueError):
        parse_regions([("chr1", "abracadabra", 5)])
        parse_regions([("ch1", 1, 2, "chr1:1-2", "puppies")])  # puppies are not allowed
        parse_regions([("chr1", 3, "abracadabra")])
        parse_regions([("chr1", 5, None)])
        parse_regions([("chr1", 5, None)], chromsizes={"chr2": 40})
        parse_regions([("chr1", 5, 0)])
