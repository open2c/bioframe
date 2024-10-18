import numpy as np
import pandas as pd
import pytest

import bioframe
from bioframe.core import specs


def test_get_default_colnames():
    assert specs._get_default_colnames() == ("chrom", "start", "end")


def test_update_default_colnames():
    new_names = ("C", "chromStart", "chromStop")
    specs.update_default_colnames(new_names)
    assert specs._get_default_colnames() == new_names

    # test that with updated default column names, bioframe.ops recognizes df1
    df1 = pd.DataFrame(
        [["chr1", 1, 5], ["chr1", 3, 8], ["chr1", 8, 10], ["chr1", 12, 14]],
        columns=list(new_names),
    )
    df1_chromsizes = {"chr1": 100, "chrX": 100}

    df1_complement = pd.DataFrame(
        [
            ["chr1", 0, 1, "chr1"],
            ["chr1", 10, 12, "chr1"],
            ["chr1", 14, 100, "chr1"],
            ["chrX", 0, 100, "chrX"],
        ],
        columns=[*list(new_names), "view_region"],
    )

    pd.testing.assert_frame_equal(
        bioframe.complement(df1, view_df=df1_chromsizes), df1_complement
    )

    # cannot update with just two colujmns
    with pytest.raises(ValueError):
        specs.update_default_colnames(("chromosome", "position"))

    # extra stuff is not allowed
    with pytest.raises(ValueError):
        specs.update_default_colnames(["chromosome", "start", "end", "extrasuff"])

    # reset to default
    specs.update_default_colnames(("chrom", "start", "end"))


def test_verify_columns():
    new_names = ("C", "chromStart", "chromStop")
    df1 = pd.DataFrame(
        [["chr1", 1, 5], ["chr1", 3, 8], ["chr1", 8, 10], ["chr1", 12, 14]],
        columns=list(new_names),
    )

    with pytest.raises(ValueError):
        specs._verify_columns(df1, specs._get_default_colnames())

    assert specs._verify_columns(
        df1,
        new_names,
        return_as_bool=True,
    )

    # no repeated column names
    with pytest.raises(ValueError):
        specs._verify_columns(df1, ["chromStart", "chromStart"], unique_cols=True)


def test_verify_column_dtypes():
    new_names = ("C", "chromStart", "chromStop")
    df1 = pd.DataFrame(
        [["chr1", 1, 5], ["chr1", 3, 8], ["chr1", 8, 10], ["chr1", 12, 14]],
        columns=list(new_names),
    )

    with pytest.raises(ValueError):
        specs._verify_column_dtypes(df1, specs._get_default_colnames())

    assert specs._verify_column_dtypes(df1, new_names, return_as_bool=True)

    df1["chromStart"] = df1["chromStart"].astype(float)
    assert not specs._verify_column_dtypes(df1, new_names, return_as_bool=True)

    df1["chromStart"] = df1["chromStart"].astype(pd.Int64Dtype())
    assert specs._verify_column_dtypes(df1, new_names, return_as_bool=True)

    df1["C"] = df1["C"].str.replace("chr", "").astype(np.int64)
    assert not specs._verify_column_dtypes(df1, new_names, return_as_bool=True)


def test_is_chrom_dtype():
    assert specs.is_chrom_dtype(str)
    fruit = pd.CategoricalDtype(
        categories=["oranges", "grapefruit", "apples"], ordered=True
    )
    assert specs.is_chrom_dtype(fruit)
    assert not specs.is_chrom_dtype(int)
    assert not specs.is_chrom_dtype(float)
