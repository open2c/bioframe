from io import StringIO

import pandas as pd
import numpy as np
import pytest

from bioframe.core.checks import *
from bioframe.ops import sort_bedframe


def test_is_cataloged():
    ### chr2q is not in view
    view_df = pd.DataFrame(
        [
            ["chr1", 0, 12, "chr1p"],
            ["chr1", 13, 26, "chr1q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    df = pd.DataFrame(
        [
            ["chr1", 0, 12, "chr1p"],
            ["chr1", 5, 15, "chr1p"],
            ["chr2", 13, 26, "chr2q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    assert not is_cataloged(df, view_df)

    ### chr1q is in view, df_view_col and view_name_col have funny labels.
    view_df = pd.DataFrame(
        [
            ["chr1", 0, 12, "chr1p"],
            ["chr1", 13, 26, "chr1q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "funny_name"],
    )
    df = pd.DataFrame(
        [
            ["chr1", 0, 12, "chr1p"],
            ["chr2", 13, 26, "chr1q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "funny_view_region"],
    )
    assert is_cataloged(
        df, view_df, df_view_col="funny_view_region", view_name_col="funny_name"
    )


def test_is_contained():

    view_df = pd.DataFrame(
        [
            ["chr1", 0, 12, "chr1p"],
            ["chr1", 13, 26, "chr1q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "name"],
    )

    ### not contained because chr2q is not cataloged
    df = pd.DataFrame(
        [
            ["chr1", 0, 12, "chr1p"],
            ["chr2", 13, 26, "chr2q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    assert not is_contained(df, view_df)

    ### not contained because first two intervals fall outside the view regions
    df = pd.DataFrame(
        [
            ["chr1", 14, 15, "chr1p"],
            ["chr1", -1, 1, "chr1q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    assert not is_contained(df, view_df)

    ### is contained
    df = pd.DataFrame(
        [
            ["chr1", 12, 12, "chr1p"],
            ["chr1", 13, 14, "chr1q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    assert is_contained(df, view_df)


def test_is_overlapping():
    ### interval on chr1 overlaps
    d = """chrom  start  end
         0  chr1      3    6
         1  chr1     5   10
         2  chr2    5  10"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    assert is_overlapping(df) is True

    ### adjacent intervals do not overlap
    d = """chrom  start  end
         0  chr1    3     6
         1  chr1    6    10
         2  chr2    5    10"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    assert is_overlapping(df) is False


def test_is_covering():
    ### test is_covering where an interval from df completely overlaps
    ### two different regions from view
    df1 = pd.DataFrame(
        [
            ["chr1", -5, 25],
        ],
        columns=["chrom", "start", "end"],
    )
    chromsizes = [("chr1", 0, 9, "chr1p"), ("chr1", 11, 20, "chr1q")]
    assert is_covering(df1, chromsizes) is True

    ### test is_covering where two intervals from df overlap
    ### two different regions from view
    df1 = pd.DataFrame(
        [
            ["chr1", -5, 10],
            ["chr1", 11, 12],
            ["chr1", 12, 20],
        ],
        columns=["chrom", "start", "end"],
    )
    chromsizes = [("chr1", 0, 9, "chr1p"), ("chr1", 11, 20, "chr1q")]
    assert is_covering(df1, chromsizes) is True

    ### test is_covering where two intervals from df overlap
    ### two different regions from view
    df1 = pd.DataFrame(
        [
            ["chr1", -5, 10, "chr1q"],
            ["chr1", 11, 12, "chr1q"],
            ["chr1", 12, 20, "chr1q"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    chromsizes = [("chr1", 0, 9, "chr1p"), ("chr1", 11, 20, "chr1q")]
    assert is_covering(df1, chromsizes) is True


def test_is_tiling():
    ### view region chr1p is tiled by one interval, chr1q is tiled by two
    df1 = pd.DataFrame(
        [
            ["chr1", 0, 9, "chr1p"],
            ["chr1", 11, 12, "chr1q"],
            ["chr1", 12, 20, "chr1q"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    chromsizes = [("chr1", 0, 9, "chr1p"), ("chr1", 11, 20, "chr1q")]
    assert is_tiling(df1, chromsizes) is True

    ### not contained, since (chr1,0,9) is associated with chr1q
    df1 = pd.DataFrame(
        [
            ["chr1", 0, 9, "chr1q"],
            ["chr1", 11, 12, "chr1q"],
            ["chr1", 12, 20, "chr1q"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    chromsizes = [("chr1", 0, 9, "chr1p"), ("chr1", 11, 20, "chr1q")]
    assert is_tiling(df1, chromsizes) is False

    ### not contained, contains overlaps
    df1 = pd.DataFrame(
        [
            ["chr1", 0, 9, "chr1p"],
            ["chr1", 11, 13, "chr1q"],
            ["chr1", 12, 20, "chr1q"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    chromsizes = [("chr1", 0, 9, "chr1p"), ("chr1", 11, 20, "chr1q")]
    assert is_tiling(df1, chromsizes) is False

    ### not covering
    df1 = pd.DataFrame(
        [
            ["chr1", 11, 12, "chr1q"],
            ["chr1", 12, 20, "chr1q"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    chromsizes = [("chr1", 0, 9, "chr1p"), ("chr1", 11, 20, "chr1q")]
    assert is_tiling(df1, chromsizes) is False


def test_is_bedframe():
    ##missing a column
    df1 = pd.DataFrame(
        [
            ["chr1", 11],
            ["chr1", 12],
        ],
        columns=["chrom", "start"],
    )
    assert is_bedframe(df1) is False

    ### end column has invalid dtype
    df1 = pd.DataFrame(
        [
            ["chr1", 10, "20"],
            ["chr1", 10, "12"],
        ],
        columns=["chrom", "start", "end"],
    )
    assert is_bedframe(df1) is False

    ### second interval start > ends.
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20],
            ["chr1", 15, 10],
        ],
        columns=["chrom", "start", "end"],
    )
    assert is_bedframe(df1) is False

    ### third interval has a null in one column
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20, "first"],
            ["chr1", 10, 15, "second"],
            ["chr1", pd.NA, 15, "third"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    # should raise  a TypeError if the second column is an object
    with pytest.raises(TypeError):
        is_bedframe(df1, raise_errors=True)
    # should raise  a ValueError after recasting to pd.Int64Dtype
    df1 = df1.astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})
    with pytest.raises(ValueError):
        is_bedframe(df1, raise_errors=True)

    ### first interval is completely NA
    df1 = pd.DataFrame(
        [
            [pd.NA, pd.NA, pd.NA, "first"],
            ["chr1", 10, 15, "second"],
            ["chr1", 10, 15, "third"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    df1 = df1.astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})
    assert is_bedframe(df1) is True


def test_is_viewframe():
    # not a bedframe
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20, "chr1p"],
            ["chr1", 15, 10, "chr1q"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    assert is_viewframe(df1) is False

    # no column for region name
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20],
            ["chr1", 30, 40],
        ],
        columns=["chrom", "start", "end"],
    )
    assert is_viewframe(df1) is False

    # contains null values
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20, "chr1p"],
            ["chr1", pd.NA, np.nan, "chr1q"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    assert is_viewframe(df1) is False

    # overlapping intervals
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20, "chr1p"],
            ["chr1", 15, 25, "chr1q"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    assert is_viewframe(df1) is False

    # valid view
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20, "chr1p"],
            ["chr1", 20, 25, "chr1q"],
            ["chr2", 20, 25, "chrTEST_2p"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    assert is_viewframe(df1) is True


def test_is_sorted():

    view_df = pd.DataFrame(
        [
            ["chrX", 1, 8, "oranges"],
            ["chrX", 8, 20, "grapefruit"],
            ["chr1", 0, 10, "apples"],
        ],
        columns=["chrom", "start", "end", "fruit"],
    )
    df_view_cat = pd.CategoricalDtype(
        categories=["oranges", "grapefruit", "apples"], ordered=True
    )
    view_df = view_df.astype({"fruit": df_view_cat})

    assert is_sorted(
        view_df, view_df=view_df, view_name_col="fruit", df_view_col="fruit"
    )

    df = pd.DataFrame(
        [
            ["chr1", 0, 10, "+"],
            ["chrX", 5, 10, "+"],
            ["chrX", 0, 5, "+"],
            ["chr2", 5, 10, "+"],
        ],
        columns=["chrom", "start", "end", "strand"],
    )

    assert not is_sorted(df)

    bfs = sort_bedframe(
        df, view_df=view_df, view_name_col="fruit")
    
    assert is_sorted(bfs, view_df=view_df, view_name_col="fruit")

    # view_df specifies a different ordering, so should not be sorted
    assert not is_sorted(bfs)
