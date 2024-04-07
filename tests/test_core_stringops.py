import numpy as np
import pytest

from bioframe.core import stringops
from bioframe.core.stringops import parse_region


def test_to_ucsc_string():
    assert stringops.to_ucsc_string(("chr21", 1, 4)) == "chr21:1-4"


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


def test_parse_region_string():
    assert stringops.parse_region_string("6:1kb-2kb") == ("6", 1000, 2000)
    assert stringops.parse_region_string("6:1,000-2,000") == ("6", 1000, 2000)
    assert stringops.parse_region_string("c6:1000-2000") == ("c6", 1000, 2000)


def test_is_complete_ucsc_string():
    assert stringops.is_complete_ucsc_string("chrX:1M-2M")
    assert not stringops.is_complete_ucsc_string("chrX")
    assert not stringops.is_complete_ucsc_string("1M-2M")
    assert not stringops.is_complete_ucsc_string(1000)
    assert not stringops.is_complete_ucsc_string(np.array([100, 200]))
    assert not stringops.is_complete_ucsc_string(np.array(["chr1:100-200"]))
