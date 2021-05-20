import pandas as pd
import numpy as np
from io import StringIO
import bioframe.core
import bioframe
import os.path as op
import pytest

testdir = op.realpath(op.dirname(__file__))


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
    assert not bioframe.core.is_cataloged(df, view_df)

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
    assert bioframe.core.is_cataloged(
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
    assert not bioframe.core.is_contained(df, view_df)

    ### not contained because first two intervals fall outside the view regions
    df = pd.DataFrame(
        [
            ["chr1", 14, 15, "chr1p"],
            ["chr1", -1, 1, "chr1q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    assert not bioframe.core.is_contained(df, view_df)

    ### is contained
    df = pd.DataFrame(
        [
            ["chr1", 12, 12, "chr1p"],
            ["chr1", 13, 14, "chr1q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    assert bioframe.core.is_contained(df, view_df)


def test_is_overlapping():
    ### interval on chr1 overlaps
    d = """chrom  start  end
         0  chr1      3    6
         1  chr1     5   10
         2  chr2    5  10"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    assert bioframe.core.is_overlapping(df) is True

    ### adjacent intervals do not overlap
    d = """chrom  start  end
         0  chr1    3     6
         1  chr1    6    10
         2  chr2    5    10"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    assert bioframe.core.is_overlapping(df) is False


def test_is_covering():
    ### test is_covering where an interval from df completely overlaps
    ### two different regions from view
    df1 = pd.DataFrame(
        [
            ["chr1", -5, 25],
        ],
        columns=["chrom", "start", "end"],
    )
    chromsizes = [ ('chr1',0,9,"chr1p"), ('chr1',11,20, "chr1q")]
    assert bioframe.core.is_covering(df1, chromsizes) is True

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
    chromsizes = [ ('chr1',0,9,"chr1p"), ('chr1',11,20, "chr1q")]
    assert bioframe.core.is_covering(df1, chromsizes) is True

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
    chromsizes = [ ('chr1',0,9,"chr1p"), ('chr1',11,20, "chr1q")]
    assert bioframe.core.is_covering(df1, chromsizes) is True


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
    chromsizes = [ ('chr1',0,9,"chr1p"), ('chr1',11,20, "chr1q")]
    assert bioframe.core.is_tiling(df1, chromsizes) is True

    ### not contained, since (chr1,0,9) is associated with chr1q
    df1 = pd.DataFrame(
        [
            ["chr1", 0, 9, "chr1q"],
            ["chr1", 11, 12, "chr1q"],
            ["chr1", 12, 20, "chr1q"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    chromsizes = [ ('chr1',0,9,"chr1p"), ('chr1',11,20, "chr1q")]
    assert bioframe.core.is_tiling(df1, chromsizes) is False

    ### not contained, contains overlaps
    df1 = pd.DataFrame(
        [
            ["chr1", 0, 9, "chr1p"],
            ["chr1", 11, 13, "chr1q"],
            ["chr1", 12, 20, "chr1q"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    chromsizes = [ ('chr1', 0, 9,"chr1p"), ('chr1', 11, 20, "chr1q")]
    assert bioframe.core.is_tiling(df1, chromsizes) is False

    ### not covering
    df1 = pd.DataFrame(
        [
            ["chr1", 11, 12, "chr1q"],
            ["chr1", 12, 20, "chr1q"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    chromsizes = [ ('chr1',0,9,"chr1p"), ('chr1',11,20, "chr1q")]
    assert bioframe.core.is_tiling(df1, chromsizes) is False


def test_is_bedframe():
    ##missing a column
    df1 = pd.DataFrame(
        [
            ["chr1", 11],
            ["chr1", 12],
        ],
        columns=["chrom", "start"],
    )
    assert bioframe.core.is_bedframe(df1) is False

    ### end column has invalid dtype
    df1 = pd.DataFrame(
        [
            ["chr1", 10, "20"],
            ["chr1", 10, "12"],
        ],
        columns=["chrom", "start", "end"],
    )
    assert bioframe.core.is_bedframe(df1) is False

    ### second interval start > ends.
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20],
            ["chr1", 15, 10],
        ],
        columns=["chrom", "start", "end"],
    )
    assert bioframe.core.is_bedframe(df1) is False


def test_sanitize_bedframe():
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20],
            ["chr1", 10, 20],
            ["chr1", 15, np.nan],
            ["chr1", pd.NA, 25],
        ],
        columns=["chrom", "start", "end"],
    )

    # drop rows with null values
    sanitized_df1 = pd.DataFrame(
        [["chr1", 10, 20], ["chr1", 10, 20]], columns=["chrom", "start", "end"]
    )
    sanitized_df1 = sanitized_df1.astype(
        {"chrom": str, "start": pd.Int64Dtype(), "end": pd.Int64Dtype()}
    )
    pd.testing.assert_frame_equal(
        sanitized_df1, bioframe.core.sanitize_bedframe(df1, drop_null=True)
    )

    # keep rows with null, but recast
    sanitized_df1 = df1.astype(
        {"chrom": str, "start": pd.Int64Dtype(), "end": pd.Int64Dtype()}
    )
    pd.testing.assert_frame_equal(sanitized_df1, bioframe.core.sanitize_bedframe(df1))

    # flip intervals as well as drop NA
    df1 = pd.DataFrame(
        [
            ["chr1", 20, 10],
            ["chr1", pd.NA, 25],
        ],
        columns=["chrom", "start", "end"],
    )
    sanitized_df1 = pd.DataFrame([["chr1", 10, 20]], columns=["chrom", "start", "end"])
    sanitized_df1 = sanitized_df1.astype(
        {"chrom": str, "start": pd.Int64Dtype(), "end": pd.Int64Dtype()}
    )
    pd.testing.assert_frame_equal(
        sanitized_df1,
        bioframe.core.sanitize_bedframe(
            df1, start_exceed_end_action="fLiP", drop_null=True
        ),
    )

    # flip intervals as well as drop NA
    df1 = pd.DataFrame(
        [
            ["chr1", 20, 10],
            ["chr1", pd.NA, 25],
        ],
        columns=["chrom", "start", "end"],
    )
    sanitized_df1 = pd.DataFrame([["chr1", 10, 20]], columns=["chrom", "start", "end"])
    sanitized_df1 = sanitized_df1.astype(
        {"chrom": str, "start": pd.Int64Dtype(), "end": pd.Int64Dtype()}
    )
    assert bioframe.core.sanitize_bedframe(
        df1, start_exceed_end_action="drop", drop_null=True
    ).empty


def test_is_viewframe():
    # not a bedframe
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20, "chr1p"],
            ["chr1", 15, 10, "chr1q"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    assert bioframe.core.is_viewframe(df1) is False

    # no column for region name
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20],
            ["chr1", 30, 40],
        ],
        columns=["chrom", "start", "end"],
    )
    assert bioframe.core.is_viewframe(df1) is False

    # contains null values
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20, "chr1p"],
            ["chr1", pd.NA, np.nan, "chr1q"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    assert bioframe.core.is_viewframe(df1) is False

    # overlapping intervals
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20, "chr1p"],
            ["chr1", 15, 25, "chr1q"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    assert bioframe.core.is_viewframe(df1) is False

    # valid view
    df1 = pd.DataFrame(
        [
            ["chr1", 10, 20, "chr1p"],
            ["chr1", 20, 25, "chr1q"],
            ["chr2", 20, 25, "chrTEST_2p"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    assert bioframe.core.is_viewframe(df1) is True


def test_make_viewframe():

    # test dict input
    d = """            chrom  start  end        name
    0    chrTESTX      0   10    chrTESTX
    1  chrTESTX_p      0   12  chrTESTX_p"""
    view_df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(
        view_df.copy(),
        bioframe.core.make_viewframe({"chrTESTX": 10, "chrTESTX_p": 12}),
    )

    # test list input
    region_list = [("chrTESTX", 0, 10), ("chrTESTX_p", 0, 12)]
    pd.testing.assert_frame_equal(
        view_df.copy(),
        bioframe.core.make_viewframe(region_list),
    )

    # test pd.Series input
    chromsizes = pd.Series(data=[5, 8], index=["chrTESTXq", "chrTEST_2p"])
    d = """      chrom  start  end        name
    0  chrTESTXq      0    5   chrTESTXq
    1   chrTEST_2p      0    8  chrTEST_2p"""
    view_df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(
        view_df.copy(), bioframe.core.make_viewframe(chromsizes)
    )

    d = """          chrom   start   end name
    0   chrTESTXq   0   5   chrTESTXq:0-5
    1   chrTEST_2p  0   8   chrTEST_2p:0-8"""
    view_df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(
        view_df.copy(), bioframe.core.make_viewframe(chromsizes, view_names_as_UCSC=True)
    )

    # test pd.DataFrame input
    pd.testing.assert_frame_equal(view_df.copy(), bioframe.core.make_viewframe(view_df))


def test_assign_view():

    ## default assignment case
    view_df = pd.DataFrame(
        [
            ["chr11", 1, 8, "chr11p"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    df = pd.DataFrame(
        [
            ["chr11", 0, 10, "+"],
        ],
        columns=["chrom", "start", "end", "strand"],
    )
    df_assigned = pd.DataFrame(
        [
            ["chr11", 0, 10, "+", "chr11p"],
        ],
        columns=["chrom", "start", "end", "strand", "view_region"],
    )
    df_assigned = df_assigned.astype(
        {"chrom": str, "start": pd.Int64Dtype(), "end": pd.Int64Dtype()}
    )
    pd.testing.assert_frame_equal(df_assigned, bioframe.core.assign_view(df, view_df))

    # assignment with funny view_name_col and an interval on chr2 not cataloged in the view_df
    view_df = pd.DataFrame(
        [
            ["chrX", 1, 8, "oranges"],
            ["chrX", 8, 20, "grapefruit"],
            ["chr1", 0, 10, "apples"],
        ],
        columns=["chrom", "start", "end", "fruit"],
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

    df_assigned = pd.DataFrame(
        [
            ["chr1", 0, 10, "+", "apples"],
            ["chrX", 5, 10, "+", "oranges"],
            ["chrX", 0, 5, "+", "oranges"],
        ],
        columns=["chrom", "start", "end", "strand", "funny_view_region"],
    )
    df_assigned = df_assigned.astype(
        {"chrom": str, "start": pd.Int64Dtype(), "end": pd.Int64Dtype()}
    )

    pd.testing.assert_frame_equal(
        df_assigned,
        bioframe.core.assign_view(
            df,
            view_df,
            view_name_col="fruit",
            df_view_col="funny_view_region",
            drop_unassigned=True,
        ),
    )

    ### keep the interval with NA as its region if drop_unassigned is False
    df_assigned = pd.DataFrame(
        [
            ["chr1", 0, 10, "+", "apples"],
            ["chrX", 5, 10, "+", "oranges"],
            ["chrX", 0, 5, "+", "oranges"],
            ["chr2", 5, 10, "+", pd.NA],
        ],
        columns=["chrom", "start", "end", "strand", "funny_view_region"],
    )
    df_assigned = df_assigned.astype(
        {"chrom": str, "start": pd.Int64Dtype(), "end": pd.Int64Dtype()}
    )

    pd.testing.assert_frame_equal(
        df_assigned,
        bioframe.core.assign_view(
            df,
            view_df,
            view_name_col="fruit",
            df_view_col="funny_view_region",
            drop_unassigned=False,
        ),
    )


def test_sort_bedframe():

    view_df = pd.DataFrame(
        [
            ["chrX", 1, 8, "oranges"],
            ["chrX", 8, 20, "grapefruit"],
            ["chr1", 0, 10, "apples"],
        ],
        columns=["chrom", "start", "end", "fruit"],
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

    ### sorting just by chrom,start,end
    df_sorted = pd.DataFrame(
        [
            ["chr1", 0, 10, "+"],
            ["chr2", 5, 10, "+"],
            ["chrX", 0, 5, "+"],
            ["chrX", 5, 10, "+"],
        ],
        columns=["chrom", "start", "end", "strand"],
    )

    pd.testing.assert_frame_equal(df_sorted, bioframe.core.sort_bedframe(df))

    ### drops unassigned chromosomes when infer_assignment is true
    df_sorted = pd.DataFrame(
        [
            ["chrX", 0, 5, "+", "oranges"],
            ["chrX", 5, 10, "+", "oranges"],
            ["chr1", 0, 10, "+", "apples"],
            ["chr2", 5, 10, "+", pd.NA],

        ],
        columns=["chrom", "start", "end", "strand", "view_region"],
    )
    df_view_cat = pd.CategoricalDtype(
        categories=["oranges", "grapefruit", "apples"], ordered=True
    )
    df_sorted = df_sorted.astype(
        {
            "chrom": str,
            "start": pd.Int64Dtype(),
            "end": pd.Int64Dtype(),
            "view_region": df_view_cat,
        }
    )

    pd.testing.assert_frame_equal(
        df_sorted, bioframe.core.sort_bedframe(df, view_df, view_name_col="fruit")
    )

    ### chr2 is not in the view, so if infer_assignment=False this should raise a ValueError
    assert pytest.raises(
        ValueError,
        bioframe.core.sort_bedframe,
        df,
        view_df,
        view_name_col="fruit",
        infer_assignment=False,
    )


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

    assert bioframe.core.is_sorted(
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

    assert not bioframe.core.is_sorted(df)

    bfs = bioframe.core.sort_bedframe(
        df, view_df=view_df, view_name_col="fruit", infer_assignment=True
    )

    assert bioframe.core.is_sorted(bfs, view_df=view_df, view_name_col="fruit")

    # view_df specifies a different ordering, so should not be sorted
    assert not bioframe.core.is_sorted(bfs)
