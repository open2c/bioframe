from io import StringIO

import pandas as pd
import numpy as np
import pytest

import bioframe

import bioframe.core.checks as checks

# import pyranges as pr

# def bioframe_to_pyranges(df):
#     pydf = df.copy()
#     pydf.rename(
#         {"chrom": "Chromosome", "start": "Start", "end": "End"},
#         axis="columns",
#         inplace=True,
#     )
#     return pr.PyRanges(pydf)


# def pyranges_to_bioframe(pydf):
#     df = pydf.df
#     df.rename(
#         {"Chromosome": "chrom", "Start": "start", "End": "end", "Count": "n_intervals"},
#         axis="columns",
#         inplace=True,
#     )
#     return df


# def pyranges_overlap_to_bioframe(pydf):
#     ## convert the df output by pyranges join into a bioframe-compatible format
#     df = pydf.df.copy()
#     df.rename(
#         {
#             "Chromosome": "chrom_1",
#             "Start": "start_1",
#             "End": "end_1",
#             "Start_b": "start_2",
#             "End_b": "end_2",
#         },
#         axis="columns",
#         inplace=True,
#     )
#     df["chrom_1"] = df["chrom_1"].values.astype("object")  # to remove categories
#     df["chrom_2"] = df["chrom_1"].values
#     return df


chroms = ["chr12", "chrX"]


def mock_bioframe(num_entries=100):
    pos = np.random.randint(1, 1e7, size=(num_entries, 2))
    df = pd.DataFrame()
    df["chrom"] = np.random.choice(chroms, num_entries)
    df["start"] = np.min(pos, axis=1)
    df["end"] = np.max(pos, axis=1)
    df.sort_values(["chrom", "start"], inplace=True)
    return df


############# tests #####################
def test_select():
    df1 = pd.DataFrame(
        [["chrX", 3, 8], ["chr1", 4, 5], ["chrX", 1, 5]],
        columns=["chrom", "start", "end"],
    )

    region1 = "chr1:4-10"
    df_result = pd.DataFrame([["chr1", 4, 5]], columns=["chrom", "start", "end"])
    pd.testing.assert_frame_equal(
        df_result, bioframe.select(df1, region1).reset_index(drop=True)
    )

    region1 = "chrX"
    df_result = pd.DataFrame(
        [["chrX", 3, 8], ["chrX", 1, 5]], columns=["chrom", "start", "end"]
    )
    pd.testing.assert_frame_equal(
        df_result, bioframe.select(df1, region1).reset_index(drop=True)
    )

    region1 = "chrX:4-6"
    df_result = pd.DataFrame(
        [["chrX", 3, 8], ["chrX", 1, 5]], columns=["chrom", "start", "end"]
    )
    pd.testing.assert_frame_equal(
        df_result, bioframe.select(df1, region1).reset_index(drop=True)
    )

    ### select with non-standard column names
    region1 = "chrX:4-6"
    new_names = ["chr", "chrstart", "chrend"]
    df1 = pd.DataFrame(
        [["chrX", 3, 8], ["chr1", 4, 5], ["chrX", 1, 5]],
        columns=new_names,
    )
    df_result = pd.DataFrame(
        [["chrX", 3, 8], ["chrX", 1, 5]],
        columns=new_names,
    )
    pd.testing.assert_frame_equal(
        df_result, bioframe.select(df1, region1, cols=new_names).reset_index(drop=True)
    )
    region1 = "chrX"
    pd.testing.assert_frame_equal(
        df_result, bioframe.select(df1, region1, cols=new_names).reset_index(drop=True)
    )

    ### select from a DataFrame with NaNs
    colnames = ["chrom", "start", "end", "view_region"]
    df = pd.DataFrame(
        [
            ["chr1", -6, 12, "chr1p"],
            [pd.NA, pd.NA, pd.NA, "chr1q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=colnames,
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})
    df_result = pd.DataFrame(
        [["chr1", -6, 12, "chr1p"]],
        columns=colnames,
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    region1 = "chr1:0-1"
    pd.testing.assert_frame_equal(
        df_result, bioframe.select(df, region1).reset_index(drop=True)
    )


def test_trim():

    ### trim with view_df
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
            ["chr1", -6, 12, "chr1p"],
            ["chr1", 0, 12, "chr1p"],
            ["chr1", 32, 36, "chr1q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    df_trimmed = pd.DataFrame(
        [
            ["chr1", 0, 12, "chr1p"],
            ["chr1", 0, 12, "chr1p"],
            ["chr1", 26, 26, "chr1q"],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    pd.testing.assert_frame_equal(df_trimmed, bioframe.trim(df, view_df=view_df))

    ### trim with view_df interpreted from dictionary for chromsizes
    chromsizes = {"chr1": 20, "chrX_0": 5}
    df = pd.DataFrame(
        [
            ["chr1", 0, 12],
            ["chr1", 13, 26],
            ["chrX_0", 1, 8],
        ],
        columns=["chrom", "startFunky", "end"],
    )
    df_trimmed = pd.DataFrame(
        [
            ["chr1", 0, 12],
            ["chr1", 13, 20],
            ["chrX_0", 1, 5],
        ],
        columns=["chrom", "startFunky", "end"],
    )
    pd.testing.assert_frame_equal(
        df_trimmed,
        bioframe.trim(
            df,
            view_df=chromsizes,
            df_view_col="chrom",
            cols=["chrom", "startFunky", "end"],
        ),
    )

    ### trim with default limits=None and negative values
    df = pd.DataFrame(
        [
            ["chr1", -4, 12, "chr1p"],
            ["chr1", 13, 26, "chr1q"],
            ["chrX", -5, -1, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "region"],
    )
    df_trimmed = pd.DataFrame(
        [
            ["chr1", 0, 12, "chr1p"],
            ["chr1", 13, 26, "chr1q"],
            ["chrX", 0, 0, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "region"],
    )
    pd.testing.assert_frame_equal(df_trimmed, bioframe.trim(df))

    ### trim when there are NaN intervals
    df = pd.DataFrame(
        [
            ["chr1", -4, 12, "chr1p"],
            [pd.NA, pd.NA, pd.NA, "chr1q"],
            ["chrX", -5, -1, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "region"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})
    df_trimmed = pd.DataFrame(
        [
            ["chr1", 0, 12, "chr1p"],
            [pd.NA, pd.NA, pd.NA, "chr1q"],
            ["chrX", 0, 0, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "region"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})
    pd.testing.assert_frame_equal(df_trimmed, bioframe.trim(df))

    ### trim with view_df and NA intervals
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
            ["chr1", -6, 12, "chr1p"],
            ["chr1", 0, 12, "chr1p"],
            [pd.NA, pd.NA, pd.NA, pd.NA],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})
    df_trimmed = pd.DataFrame(
        [
            ["chr1", 0, 12, "chr1p"],
            ["chr1", 0, 12, "chr1p"],
            [pd.NA, pd.NA, pd.NA, pd.NA],
            ["chrX", 1, 8, "chrX_0"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})
    pd.testing.assert_frame_equal(df_trimmed, bioframe.trim(df, view_df=view_df))


def test_expand():

    d = """chrom  start  end
         0  chr1      1    5
         1  chr1     50   55
         2  chr2    100  200"""
    fake_bioframe = pd.read_csv(StringIO(d), sep=r"\s+")

    expand_bp = 10
    fake_expanded = bioframe.expand(fake_bioframe, expand_bp)
    d = """chrom  start  end
         0  chr1      -9    15
         1  chr1     40   65
         2  chr2    90  210"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, fake_expanded)

    # expand with negative pad
    expand_bp = -10
    fake_expanded = bioframe.expand(fake_bioframe, expand_bp)
    d = """chrom  start  end
         0  chr1      3    3
         1  chr1     52   52
         2  chr2    110  190"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, fake_expanded)

    expand_bp = -10
    fake_expanded = bioframe.expand(fake_bioframe, expand_bp, side="left")
    d = """chrom  start  end
         0  chr1      3    5
         1  chr1     52   55
         2  chr2    110  200"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, fake_expanded)

    # expand with multiplicative pad
    mult = 0
    fake_expanded = bioframe.expand(fake_bioframe, pad=None, scale=mult)
    d = """chrom  start  end
         0  chr1      3    3
         1  chr1     52   52
         2  chr2    150  150"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, fake_expanded)

    mult = 2.0
    fake_expanded = bioframe.expand(fake_bioframe, pad=None, scale=mult)
    d = """chrom  start  end
         0  chr1      -1    7
         1  chr1     48   58
         2  chr2    50  250"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, fake_expanded)

    # expand with NA and non-integer multiplicative pad
    d = """chrom  start  end
         0  chr1      1    5
         1  NA     NA   NA
         2  chr2    100  200"""
    fake_bioframe = pd.read_csv(StringIO(d), sep=r"\s+").astype(
        {"start": pd.Int64Dtype(), "end": pd.Int64Dtype()}
    )
    mult = 1.10
    fake_expanded = bioframe.expand(fake_bioframe, pad=None, scale=mult)
    d = """chrom  start  end
         0  chr1      1   5
         1  NA     NA   NA
         2  chr2    95  205"""
    df = pd.read_csv(StringIO(d), sep=r"\s+").astype(
        {"start": pd.Int64Dtype(), "end": pd.Int64Dtype()}
    )
    pd.testing.assert_frame_equal(df, fake_expanded)


def test_overlap():

    ### test consistency of overlap(how='inner') with pyranges.join ###
    ### note does not test overlap_start or overlap_end columns of bioframe.overlap
    df1 = mock_bioframe()
    df2 = mock_bioframe()
    assert df1.equals(df2) == False

    # p1 = bioframe_to_pyranges(df1)
    # p2 = bioframe_to_pyranges(df2)
    # pp = pyranges_overlap_to_bioframe(p1.join(p2, how=None))[
    #     ["chrom_1", "start_1", "end_1", "chrom_2", "start_2", "end_2"]
    # ]
    # bb = bioframe.overlap(df1, df2, how="inner")[
    #     ["chrom_1", "start_1", "end_1", "chrom_2", "start_2", "end_2"]
    # ]
    # pp = pp.sort_values(
    #     ["chrom_1", "start_1", "end_1", "chrom_2", "start_2", "end_2"],
    #     ignore_index=True,
    # )
    # bb = bb.sort_values(
    #     ["chrom_1", "start_1", "end_1", "chrom_2", "start_2", "end_2"],
    #     ignore_index=True,
    # )
    # pd.testing.assert_frame_equal(bb, pp, check_dtype=False, check_exact=False)
    # print("overlap elements agree")

    ### test overlap on= [] ###
    df1 = pd.DataFrame(
        [
            ["chr1", 8, 12, "+", "cat"],
            ["chr1", 8, 12, "-", "cat"],
            ["chrX", 1, 8, "+", "cat"],
        ],
        columns=["chrom1", "start", "end", "strand", "animal"],
    )

    df2 = pd.DataFrame(
        [["chr1", 6, 10, "+", "dog"], ["chrX", 7, 10, "-", "dog"]],
        columns=["chrom2", "start2", "end2", "strand", "animal"],
    )

    b = bioframe.overlap(
        df1,
        df2,
        on=["animal"],
        how="left",
        cols1=("chrom1", "start", "end"),
        cols2=("chrom2", "start2", "end2"),
        return_index=True,
        return_input=False,
    )
    assert np.sum(pd.isna(b["index_"].values)) == 3

    b = bioframe.overlap(
        df1,
        df2,
        on=["strand"],
        how="left",
        cols1=("chrom1", "start", "end"),
        cols2=("chrom2", "start2", "end2"),
        return_index=True,
        return_input=False,
    )
    assert np.sum(pd.isna(b["index_"].values)) == 2

    b = bioframe.overlap(
        df1,
        df2,
        on=None,
        how="left",
        cols1=("chrom1", "start", "end"),
        cols2=("chrom2", "start2", "end2"),
        return_index=True,
        return_input=False,
    )
    assert np.sum(pd.isna(b["index_"].values)) == 0

    ### test overlap 'left', 'outer', and 'right'
    b = bioframe.overlap(
        df1,
        df2,
        on=None,
        how="outer",
        cols1=("chrom1", "start", "end"),
        cols2=("chrom2", "start2", "end2"),
    )
    assert len(b) == 3

    b = bioframe.overlap(
        df1,
        df2,
        on=["animal"],
        how="outer",
        cols1=("chrom1", "start", "end"),
        cols2=("chrom2", "start2", "end2"),
    )
    assert len(b) == 5

    b = bioframe.overlap(
        df1,
        df2,
        on=["animal"],
        how="inner",
        cols1=("chrom1", "start", "end"),
        cols2=("chrom2", "start2", "end2"),
    )
    assert len(b) == 0

    b = bioframe.overlap(
        df1,
        df2,
        on=["animal"],
        how="right",
        cols1=("chrom1", "start", "end"),
        cols2=("chrom2", "start2", "end2"),
    )
    assert len(b) == 2

    b = bioframe.overlap(
        df1,
        df2,
        on=["animal"],
        how="left",
        cols1=("chrom1", "start", "end"),
        cols2=("chrom2", "start2", "end2"),
    )
    assert len(b) == 3

    ### test keep_order and NA handling
    df1 = pd.DataFrame(
        [
            ["chr1", 8, 12, "+"],
            [pd.NA, pd.NA, pd.NA, "-"],
            ["chrX", 1, 8, "+"],
        ],
        columns=["chrom", "start", "end", "strand"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    df2 = pd.DataFrame(
        [["chr1", 6, 10, "+"], [pd.NA, pd.NA, pd.NA, "-"], ["chrX", 7, 10, "-"]],
        columns=["chrom2", "start2", "end2", "strand"],
    ).astype({"start2": pd.Int64Dtype(), "end2": pd.Int64Dtype()})

    assert df1.equals(
        bioframe.overlap(
            df1, df2, how="left", keep_order=True, cols2=["chrom2", "start2", "end2"]
        )[["chrom", "start", "end", "strand"]]
    )
    assert ~df1.equals(
        bioframe.overlap(
            df1, df2, how="left", keep_order=False, cols2=["chrom2", "start2", "end2"]
        )[["chrom", "start", "end", "strand"]]
    )

    df1 = pd.DataFrame(
        [
            ["chr1", 8, 12, "+", pd.NA],
            [pd.NA, pd.NA, pd.NA, "-", pd.NA],
            ["chrX", 1, 8, pd.NA, pd.NA],
        ],
        columns=["chrom", "start", "end", "strand", "animal"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    df2 = pd.DataFrame(
        [["chr1", 6, 10, pd.NA, "tiger"]],
        columns=["chrom2", "start2", "end2", "strand", "animal"],
    ).astype({"start2": pd.Int64Dtype(), "end2": pd.Int64Dtype()})

    assert (
        bioframe.overlap(
            df1,
            df2,
            how="outer",
            cols2=["chrom2", "start2", "end2"],
            return_index=True,
            keep_order=False,
        ).shape
        == (3, 12)
    )

    ### result of overlap should still have bedframe-like properties
    overlap_df = bioframe.overlap(
        df1,
        df2,
        how="outer",
        cols2=["chrom2", "start2", "end2"],
        return_index=True,
        suffixes=("", ""),
    )
    assert checks.is_bedframe(
        overlap_df[df1.columns],
    )
    assert checks.is_bedframe(
        overlap_df[df2.columns], cols=["chrom2", "start2", "end2"]
    )

    overlap_df = bioframe.overlap(
        df1,
        df2,
        how="innter",
        cols2=["chrom2", "start2", "end2"],
        return_index=True,
        suffixes=("", ""),
    )
    assert checks.is_bedframe(
        overlap_df[df1.columns],
    )
    assert checks.is_bedframe(
        overlap_df[df2.columns], cols=["chrom2", "start2", "end2"]
    )

    # test keep_order incompatible if how!= 'left'
    with pytest.raises(ValueError):
        bioframe.overlap(
            df1,
            df2,
            how="outer",
            on=["animal"],
            cols2=["chrom2", "start2", "end2"],
            keep_order=True,
        )


def test_cluster():
    df1 = pd.DataFrame(
        [
            ["chr1", 1, 5],
            ["chr1", 3, 8],
            ["chr1", 8, 10],
            ["chr1", 12, 14],
        ],
        columns=["chrom", "start", "end"],
    )
    df_annotated = bioframe.cluster(df1)
    assert (
        df_annotated["cluster"].values == np.array([0, 0, 0, 1])
    ).all()  # the last interval does not overlap the first three
    df_annotated = bioframe.cluster(df1, min_dist=2)
    assert (
        df_annotated["cluster"].values == np.array([0, 0, 0, 0])
    ).all()  # all intervals part of the same cluster

    df_annotated = bioframe.cluster(df1, min_dist=None)
    assert (
        df_annotated["cluster"].values == np.array([0, 0, 1, 2])
    ).all()  # adjacent intervals not clustered

    df1.iloc[0, 0] = "chrX"
    df_annotated = bioframe.cluster(df1)
    assert (
        df_annotated["cluster"].values == np.array([2, 0, 0, 1])
    ).all()  # do not cluster intervals across chromosomes

    # test consistency with pyranges (which automatically sorts df upon creation and uses 1-based indexing for clusters)
    # assert (
    #     (bioframe_to_pyranges(df1).cluster(count=True).df["Cluster"].values - 1)
    #     == bioframe.cluster(df1.sort_values(["chrom", "start"]))["cluster"].values
    # ).all()

    # test on=[] argument
    df1 = pd.DataFrame(
        [
            ["chr1", 3, 8, "+", "cat", 5.5],
            ["chr1", 3, 8, "-", "dog", 6.5],
            ["chr1", 6, 10, "-", "cat", 6.5],
            ["chrX", 6, 10, "-", "cat", 6.5],
        ],
        columns=["chrom", "start", "end", "strand", "animal", "location"],
    )
    assert (
        bioframe.cluster(df1, on=["animal"])["cluster"].values == np.array([0, 1, 0, 2])
    ).all()
    assert (
        bioframe.cluster(df1, on=["strand"])["cluster"].values == np.array([0, 1, 1, 2])
    ).all()
    assert (
        bioframe.cluster(df1, on=["location", "animal"])["cluster"].values
        == np.array([0, 2, 1, 3])
    ).all()

    ### test cluster with NAs
    df1 = pd.DataFrame(
        [
            ["chrX", 1, 8, pd.NA, pd.NA],
            [pd.NA, pd.NA, pd.NA, "-", pd.NA],
            ["chr1", 8, 12, "+", pd.NA],
            ["chr1", 1, 8, np.nan, pd.NA],
            [pd.NA, np.nan, pd.NA, "-", pd.NA],
        ],
        columns=["chrom", "start", "end", "strand", "animal"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    assert bioframe.cluster(df1)["cluster"].max() == 3
    assert bioframe.cluster(df1, on=["strand"])["cluster"].max() == 4
    pd.testing.assert_frame_equal(df1, bioframe.cluster(df1)[df1.columns])

    assert checks.is_bedframe(
        bioframe.cluster(df1, on=["strand"]),
        cols=["chrom", "cluster_start", "cluster_end"],
    )
    assert checks.is_bedframe(
        bioframe.cluster(df1), cols=["chrom", "cluster_start", "cluster_end"]
    )
    assert checks.is_bedframe(bioframe.cluster(df1))


def test_merge():
    df1 = pd.DataFrame(
        [
            ["chr1", 1, 5],
            ["chr1", 3, 8],
            ["chr1", 8, 10],
            ["chr1", 12, 14],
        ],
        columns=["chrom", "start", "end"],
    )

    # the last interval does not overlap the first three with default min_dist=0
    assert (bioframe.merge(df1)["n_intervals"].values == np.array([3, 1])).all()

    # adjacent intervals are not clustered with min_dist=none
    assert (
        bioframe.merge(df1, min_dist=None)["n_intervals"].values == np.array([2, 1, 1])
    ).all()

    # all intervals part of one cluster
    assert (
        bioframe.merge(df1, min_dist=2)["n_intervals"].values == np.array([4])
    ).all()

    df1.iloc[0, 0] = "chrX"
    assert (
        bioframe.merge(df1, min_dist=None)["n_intervals"].values
        == np.array([1, 1, 1, 1])
    ).all()
    assert (
        bioframe.merge(df1, min_dist=0)["n_intervals"].values == np.array([2, 1, 1])
    ).all()

    # total number of intervals should equal length of original dataframe
    mock_df = mock_bioframe()
    assert np.sum(bioframe.merge(mock_df, min_dist=0)["n_intervals"].values) == len(
        mock_df
    )

    # # test consistency with pyranges
    # pd.testing.assert_frame_equal(
    #     pyranges_to_bioframe(bioframe_to_pyranges(df1).merge(count=True)),
    #     bioframe.merge(df1),
    #     check_dtype=False,
    #     check_exact=False,
    # )

    # test on=['chrom',...] argument
    df1 = pd.DataFrame(
        [
            ["chr1", 3, 8, "+", "cat", 5.5],
            ["chr1", 3, 8, "-", "dog", 6.5],
            ["chr1", 6, 10, "-", "cat", 6.5],
            ["chrX", 6, 10, "-", "cat", 6.5],
        ],
        columns=["chrom", "start", "end", "strand", "animal", "location"],
    )
    assert len(bioframe.merge(df1, on=None)) == 2
    assert len(bioframe.merge(df1, on=["strand"])) == 3
    assert len(bioframe.merge(df1, on=["strand", "location"])) == 3
    assert len(bioframe.merge(df1, on=["strand", "location", "animal"])) == 4
    d = """ chrom   start   end animal  n_intervals
        0   chr1    3   10  cat 2
        1   chr1    3   8   dog 1
        2   chrX    6   10  cat 1"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(
        df,
        bioframe.merge(df1, on=["animal"]),
        check_dtype=False,
    )

    # merge with repeated indices
    df = pd.DataFrame(
        {"chrom": ["chr1", "chr2"], "start": [100, 400], "end": [110, 410]}
    )
    df.index = [0, 0]
    pd.testing.assert_frame_equal(
        df.reset_index(drop=True), bioframe.merge(df)[["chrom", "start", "end"]]
    )

    # test merge with NAs
    df1 = pd.DataFrame(
        [
            ["chrX", 1, 8, pd.NA, pd.NA],
            [pd.NA, pd.NA, pd.NA, "-", pd.NA],
            ["chr1", 8, 12, "+", pd.NA],
            ["chr1", 1, 8, np.nan, pd.NA],
            [pd.NA, np.nan, pd.NA, "-", pd.NA],
        ],
        columns=["chrom", "start", "end", "strand", "animal"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    assert bioframe.merge(df1).shape[0] == 4
    assert bioframe.merge(df1)["start"].iloc[0] == 1
    assert bioframe.merge(df1)["end"].iloc[0] == 12
    assert bioframe.merge(df1, on=["strand"]).shape[0] == df1.shape[0]
    assert bioframe.merge(df1, on=["animal"]).shape[0] == df1.shape[0]
    assert bioframe.merge(df1, on=["animal"]).shape[1] == df1.shape[1] + 1
    assert checks.is_bedframe(bioframe.merge(df1, on=["strand", "animal"]))


def test_complement():
    ### complementing a df with no intervals in chrX by a view with chrX should return entire chrX region
    df1 = pd.DataFrame(
        [["chr1", 1, 5], ["chr1", 3, 8], ["chr1", 8, 10], ["chr1", 12, 14]],
        columns=["chrom", "start", "end"],
    )
    df1_chromsizes = {"chr1": 100, "chrX": 100}

    df1_complement = pd.DataFrame(
        [
            ["chr1", 0, 1, "chr1"],
            ["chr1", 10, 12, "chr1"],
            ["chr1", 14, 100, "chr1"],
            ["chrX", 0, 100, "chrX"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )

    pd.testing.assert_frame_equal(
        bioframe.complement(df1, view_df=df1_chromsizes), df1_complement
    )

    ### test complement with two chromosomes ###
    df1.iloc[0, 0] = "chrX"
    df1_complement = pd.DataFrame(
        [
            ["chr1", 0, 3, "chr1"],
            ["chr1", 10, 12, "chr1"],
            ["chr1", 14, 100, "chr1"],
            ["chrX", 0, 1, "chrX"],
            ["chrX", 5, 100, "chrX"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    pd.testing.assert_frame_equal(
        bioframe.complement(df1, view_df=df1_chromsizes), df1_complement
    )

    ### test complement with no view_df and a negative interval
    df1 = pd.DataFrame(
        [["chr1", -5, 5], ["chr1", 10, 20]], columns=["chrom", "start", "end"]
    )
    df1_complement = pd.DataFrame(
        [["chr1", 5, 10, "chr1"], ["chr1", 20, np.iinfo(np.int64).max, "chr1"]],
        columns=["chrom", "start", "end", "view_region"],
    )
    pd.testing.assert_frame_equal(bioframe.complement(df1), df1_complement)

    ### test complement with an overhanging interval
    df1 = pd.DataFrame(
        [["chr1", -5, 5], ["chr1", 10, 20]], columns=["chrom", "start", "end"]
    )
    chromsizes = {"chr1": 15}
    df1_complement = pd.DataFrame(
        [
            ["chr1", 5, 10, "chr1"],
        ],
        columns=["chrom", "start", "end", "view_region"],
    )
    pd.testing.assert_frame_equal(
        bioframe.complement(df1, view_df=chromsizes, view_name_col="VR"), df1_complement
    )

    ### test complement where an interval from df overlaps two different regions from view
    ### test complement with no view_df and a negative interval
    df1 = pd.DataFrame([["chr1", 5, 15]], columns=["chrom", "start", "end"])
    chromsizes = [("chr1", 0, 9, "chr1p"), ("chr1", 11, 20, "chr1q")]
    df1_complement = pd.DataFrame(
        [["chr1", 0, 5, "chr1p"], ["chr1", 15, 20, "chr1q"]],
        columns=["chrom", "start", "end", "view_region"],
    )
    pd.testing.assert_frame_equal(bioframe.complement(df1, chromsizes), df1_complement)

    ### test complement with NAs
    df1 = pd.DataFrame(
        [[pd.NA, pd.NA, pd.NA], ["chr1", 5, 15], [pd.NA, pd.NA, pd.NA]],
        columns=["chrom", "start", "end"],
    ).astype(
        {
            "start": pd.Int64Dtype(),
            "end": pd.Int64Dtype(),
        }
    )

    pd.testing.assert_frame_equal(bioframe.complement(df1, chromsizes), df1_complement)

    with pytest.raises(ValueError):  # no NAs allowed in chromsizes
        bioframe.complement(
            df1, [("chr1", pd.NA, 9, "chr1p"), ("chr1", 11, 20, "chr1q")]
        )

    assert checks.is_bedframe(bioframe.complement(df1, chromsizes))


def test_closest():
    df1 = pd.DataFrame(
        [
            ["chr1", 1, 5],
        ],
        columns=["chrom", "start", "end"],
    )

    df2 = pd.DataFrame(
        [["chr1", 4, 8], ["chr1", 10, 11]], columns=["chrom", "start", "end"]
    )

    ### closest(df1,df2,k=1) ###
    d = """chrom  start  end chrom_  start_  end_  distance
        0    chr1        1      5    chr1        4      8         0"""
    df = pd.read_csv(StringIO(d), sep=r"\s+").astype(
        {
            "start_": pd.Int64Dtype(),
            "end_": pd.Int64Dtype(),
            "distance": pd.Int64Dtype(),
        }
    )
    pd.testing.assert_frame_equal(df, bioframe.closest(df1, df2, k=1))

    ### closest(df1,df2, ignore_overlaps=True)) ###
    d = """chrom_1 start_1 end_1   chrom_2 start_2 end_2   distance
        0   chr1    1   5   chr1    10  11  5"""
    df = pd.read_csv(StringIO(d), sep=r"\s+").astype(
        {
            "start_2": pd.Int64Dtype(),
            "end_2": pd.Int64Dtype(),
            "distance": pd.Int64Dtype(),
        }
    )
    pd.testing.assert_frame_equal(
        df, bioframe.closest(df1, df2, suffixes=("_1", "_2"), ignore_overlaps=True)
    )

    ### closest(df1,df2,k=2) ###
    d = """chrom_1 start_1 end_1   chrom_2 start_2 end_2   distance
            0   chr1    1   5   chr1    4   8   0
            1   chr1    1   5   chr1    10  11  5"""
    df = pd.read_csv(StringIO(d), sep=r"\s+").astype(
        {
            "start_2": pd.Int64Dtype(),
            "end_2": pd.Int64Dtype(),
            "distance": pd.Int64Dtype(),
        }
    )
    pd.testing.assert_frame_equal(
        df, bioframe.closest(df1, df2, suffixes=("_1", "_2"), k=2)
    )

    ### closest(df2,df1) ###
    d = """chrom_1  start_1 end_1   chrom_2 start_2 end_2   distance
            0   chr1    4   8   chr1    1   5   0
            1   chr1    10  11  chr1    1   5   5 """
    df = pd.read_csv(StringIO(d), sep=r"\s+").astype(
        {
            "start_2": pd.Int64Dtype(),
            "end_2": pd.Int64Dtype(),
            "distance": pd.Int64Dtype(),
        }
    )
    pd.testing.assert_frame_equal(df, bioframe.closest(df2, df1, suffixes=("_1", "_2")))

    ### change first interval to new chrom ###
    df2.iloc[0, 0] = "chrA"
    d = """chrom start   end     chrom_ start_ end_  distance
              0   chr1    1   5   chr1    10  11  5"""
    df = pd.read_csv(StringIO(d), sep=r"\s+").astype(
        {
            "start_": pd.Int64Dtype(),
            "end_": pd.Int64Dtype(),
            "distance": pd.Int64Dtype(),
        }
    )
    pd.testing.assert_frame_equal(df, bioframe.closest(df1, df2, k=1))

    ### test other return arguments ###
    df2.iloc[0, 0] = "chr1"
    d = """
        index index_ have_overlap    overlap_start   overlap_end distance
        0   0   0   True    4   5   0
        1   0   1   False   <NA>    <NA>    5
        """
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(
        df,
        bioframe.closest(
            df1,
            df2,
            k=2,
            return_overlap=True,
            return_index=True,
            return_input=False,
            return_distance=True,
        ),
        check_dtype=False,
    )

    # closest should ignore empty groups (e.g. from categorical chrom)
    df = pd.DataFrame(
        [
            ["chrX", 1, 8],
            ["chrX", 2, 10],
        ],
        columns=["chrom", "start", "end"],
    )
    d = """ chrom_1 start_1 end_1   chrom_2 start_2 end_2   distance
            0   chrX    1   8   chrX    2   10  0
            1   chrX    2   10  chrX    1   8   0"""
    df_closest = pd.read_csv(StringIO(d), sep=r"\s+")
    df_cat = pd.CategoricalDtype(categories=["chrX", "chr1"], ordered=True)
    df = df.astype({"chrom": df_cat})
    pd.testing.assert_frame_equal(
        df_closest,
        bioframe.closest(df, suffixes=("_1", "_2")),
        check_dtype=False,
        check_categorical=False,
    )

    # closest should ignore null rows: code will need to be modified
    # as for overlap if an on=[] option is added
    df1 = pd.DataFrame(
        [
            [pd.NA, pd.NA, pd.NA],
            ["chr1", 1, 5],
        ],
        columns=["chrom", "start", "end"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    df2 = pd.DataFrame(
        [
            [pd.NA, pd.NA, pd.NA],
            ["chr1", 4, 8],
            [pd.NA, pd.NA, pd.NA],
            ["chr1", 10, 11],
        ],
        columns=["chrom", "start", "end"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    d = """chrom_1 start_1 end_1   chrom_2 start_2 end_2   distance
        0   chr1    1   5   chr1    10  11  5"""
    df = pd.read_csv(StringIO(d), sep=r"\s+").astype(
        {
            "start_1": pd.Int64Dtype(),
            "end_1": pd.Int64Dtype(),
            "start_2": pd.Int64Dtype(),
            "end_2": pd.Int64Dtype(),
            "distance": pd.Int64Dtype(),
        }
    )
    pd.testing.assert_frame_equal(
        df, bioframe.closest(df1, df2, suffixes=("_1", "_2"), ignore_overlaps=True, k=5)
    )

    with pytest.raises(ValueError):  # inputs must be valid bedFrames
        df1.iloc[0, 0] = "chr10"
        bioframe.closest(df1, df2)


def test_coverage():

    #### coverage does not exceed length of original interval
    df1 = pd.DataFrame([["chr1", 3, 8]], columns=["chrom", "start", "end"])
    df2 = pd.DataFrame([["chr1", 2, 10]], columns=["chrom", "start", "end"])
    d = """chrom    start   end coverage
         0  chr1    3   8   5"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    print(df)
    print(bioframe.coverage(df1, df2))
    pd.testing.assert_frame_equal(df, bioframe.coverage(df1, df2))

    ### coverage of interval on different chrom returns zero for coverage and n_overlaps
    df1 = pd.DataFrame([["chr1", 3, 8]], columns=["chrom", "start", "end"])
    df2 = pd.DataFrame([["chrX", 3, 8]], columns=["chrom", "start", "end"])
    d = """chrom    start   end coverage
        0  chr1      3       8     0   """
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, bioframe.coverage(df1, df2))

    ### when a second overlap starts within the first
    df1 = pd.DataFrame([["chr1", 3, 8]], columns=["chrom", "start", "end"])
    df2 = pd.DataFrame(
        [["chr1", 3, 6], ["chr1", 5, 8]], columns=["chrom", "start", "end"]
    )

    d = """chrom    start   end coverage
         0  chr1     3       8     5"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, bioframe.coverage(df1, df2))


def test_subtract():
    ### no intervals should be left after self-subtraction
    df1 = pd.DataFrame(
        [["chrX", 3, 8], ["chr1", 4, 7], ["chrX", 1, 5]],
        columns=["chrom", "start", "end"],
    )
    assert len(bioframe.subtract(df1, df1)) == 0

    ### no intervals on chrX should remain after subtracting a longer interval
    ### interval on chr1 should be split.
    ### additional column should be propagated to children.
    df2 = pd.DataFrame(
        [
            ["chrX", 0, 18],
            ["chr1", 5, 6],
        ],
        columns=["chrom", "start", "end"],
    )
    df1["animal"] = "sea-creature"
    df_result = pd.DataFrame(
        [["chr1", 4, 5, "sea-creature"], ["chr1", 6, 7, "sea-creature"]],
        columns=["chrom", "start", "end", "animal"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})
    pd.testing.assert_frame_equal(
        df_result,
        bioframe.subtract(df1, df2)
        .sort_values(["chrom", "start", "end"])
        .reset_index(drop=True),
    )

    ### no intervals on chrX should remain after subtracting a longer interval
    df2 = pd.DataFrame(
        [["chrX", 0, 4], ["chr1", 6, 6], ["chrX", 4, 9]],
        columns=["chrom", "start", "end"],
    )

    df1["animal"] = "sea-creature"
    df_result = pd.DataFrame(
        [["chr1", 4, 6, "sea-creature"], ["chr1", 6, 7, "sea-creature"]],
        columns=["chrom", "start", "end", "animal"],
    )
    pd.testing.assert_frame_equal(
        df_result.astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()}),
        bioframe.subtract(df1, df2)
        .sort_values(["chrom", "start", "end"])
        .reset_index(drop=True),
    )

    ### subtracting dataframes funny column names

    funny_cols = ["C", "chromStart", "chromStop"]

    df1 = pd.DataFrame(
        [["chrX", 3, 8], ["chr1", 4, 7], ["chrX", 1, 5]],
        columns=funny_cols,
    )
    df1["strand"] = "+"

    assert len(bioframe.subtract(df1, df1, cols1=funny_cols, cols2=funny_cols)) == 0

    funny_cols2 = ["chr", "st", "e"]
    df2 = pd.DataFrame(
        [
            ["chrX", 0, 18],
            ["chr1", 5, 6],
        ],
        columns=funny_cols2,
    )

    df_result = pd.DataFrame(
        [["chr1", 4, 5, "+"], ["chr1", 6, 7, "+"]],
        columns=funny_cols + ["strand"],
    )
    df_result = df_result.astype(
        {funny_cols[1]: pd.Int64Dtype(), funny_cols[2]: pd.Int64Dtype()}
    )

    pd.testing.assert_frame_equal(
        df_result,
        bioframe.subtract(df1, df2, cols1=funny_cols, cols2=funny_cols2)
        .sort_values(funny_cols)
        .reset_index(drop=True),
    )

    # subtract should ignore empty groups
    df1 = pd.DataFrame(
        [
            ["chrX", 1, 8],
            ["chrX", 2, 10],
        ],
        columns=["chrom", "start", "end"],
    )
    df2 = pd.DataFrame(
        [
            ["chrX", 1, 8],
        ],
        columns=["chrom", "start", "end"],
    )
    df_cat = pd.CategoricalDtype(categories=["chrX", "chr1"], ordered=True)
    df1 = df1.astype({"chrom": df_cat})
    df_subtracted = pd.DataFrame(
        [
            ["chrX", 8, 10],
        ],
        columns=["chrom", "start", "end"],
    )

    assert bioframe.subtract(df1, df1).empty

    pd.testing.assert_frame_equal(
        df_subtracted.astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()}),
        bioframe.subtract(df1, df2),
        check_dtype=False,
        check_categorical=False,
    )

    ## test transferred from deprecated bioframe.split
    df1 = pd.DataFrame(
        [["chrX", 3, 8], ["chr1", 4, 7], ["chrX", 1, 5]],
        columns=["chrom", "start", "end"],
    )

    df2 = pd.DataFrame(
        [
            ["chrX", 4],
            ["chr1", 5],
        ],
        columns=["chrom", "pos"],
    )
    df2["start"] = df2["pos"]
    df2["end"] = df2["pos"]

    df_result = (
        pd.DataFrame(
            [
                ["chrX", 1, 4],
                ["chrX", 3, 4],
                ["chrX", 4, 5],
                ["chrX", 4, 8],
                ["chr1", 5, 7],
                ["chr1", 4, 5],
            ],
            columns=["chrom", "start", "end"],
        )
        .sort_values(["chrom", "start", "end"])
        .reset_index(drop=True)
        .astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})
    )

    pd.testing.assert_frame_equal(
        df_result,
        bioframe.subtract(df1, df2)
        .sort_values(["chrom", "start", "end"])
        .reset_index(drop=True),
    )

    # Test the case when a chromosome should not be split (now implemented with subtract)
    df1 = pd.DataFrame(
        [
            ["chrX", 3, 8],
            ["chr1", 4, 7],
        ],
        columns=["chrom", "start", "end"],
    )

    df2 = pd.DataFrame([["chrX", 4]], columns=["chrom", "pos"])
    df2["start"] = df2["pos"].values
    df2["end"] = df2["pos"].values

    df_result = (
        pd.DataFrame(
            [
                ["chrX", 3, 4],
                ["chrX", 4, 8],
                ["chr1", 4, 7],
            ],
            columns=["chrom", "start", "end"],
        )
        .sort_values(["chrom", "start", "end"])
        .reset_index(drop=True)
    )

    pd.testing.assert_frame_equal(
        df_result.astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()}),
        bioframe.subtract(df1, df2)
        .sort_values(["chrom", "start", "end"])
        .reset_index(drop=True),
    )

    # subtract should ignore null rows
    df1 = pd.DataFrame(
        [[pd.NA, pd.NA, pd.NA], ["chr1", 1, 5]],
        columns=["chrom", "start", "end"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    df2 = pd.DataFrame(
        [
            ["chrX", 1, 5],
            [pd.NA, pd.NA, pd.NA],
            ["chr1", 4, 8],
            [pd.NA, pd.NA, pd.NA],
            ["chr1", 10, 11],
        ],
        columns=["chrom", "start", "end"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    df_subtracted = pd.DataFrame(
        [
            ["chr1", 1, 4],
        ],
        columns=["chrom", "start", "end"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    pd.testing.assert_frame_equal(df_subtracted, bioframe.subtract(df1, df2))

    df1 = pd.DataFrame(
        [
            [pd.NA, pd.NA, pd.NA],
        ],
        columns=["chrom", "start", "end"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    assert len(bioframe.subtract(df1, df2)) == 0  # empty df1 but valid chroms in df2

    with pytest.raises(ValueError):  # no non-null chromosomes
        bioframe.subtract(df1, df1)

    df2 = pd.DataFrame(
        [
            [pd.NA, pd.NA, pd.NA],
            [pd.NA, pd.NA, pd.NA],
        ],
        columns=["chrom", "start", "end"],
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})
    with pytest.raises(ValueError):  # no non-null chromosomes
        bioframe.subtract(df1, df2)


def test_setdiff():

    cols1 = ["chrom1", "start", "end"]
    cols2 = ["chrom2", "start", "end"]
    df1 = pd.DataFrame(
        [
            ["chr1", 8, 12, "+", "cat"],
            ["chr1", 8, 12, "-", "cat"],
            ["chrX", 1, 8, "+", "cat"],
        ],
        columns=cols1 + ["strand", "animal"],
    )
    df2 = pd.DataFrame(
        [
            ["chrX", 7, 10, "-", "dog"],
            ["chr1", 6, 10, "-", "cat"],
            ["chr1", 6, 10, "-", "cat"],
        ],
        columns=cols2 + ["strand", "animal"],
    )

    assert (
        len(
            bioframe.setdiff(
                df1,
                df2,
                cols1=cols1,
                cols2=cols2,
                on=None,
            )
        )
        == 0
    )  # everything overlaps

    assert (
        len(
            bioframe.setdiff(
                df1,
                df2,
                cols1=cols1,
                cols2=cols2,
                on=["animal"],
            )
        )
        == 1
    )  # two overlap, one remains

    assert (
        len(
            bioframe.setdiff(
                df1,
                df2,
                cols1=cols1,
                cols2=cols2,
                on=["strand"],
            )
        )
        == 2
    )  # one overlaps, two remain

    # setdiff should ignore nan rows
    df1 = pd.concat([pd.DataFrame([pd.NA]), df1, pd.DataFrame([pd.NA])])[
        ["chrom1", "start", "end", "strand", "animal"]
    ]
    df1 = df1.astype(
        {
            "start": pd.Int64Dtype(),
            "end": pd.Int64Dtype(),
        }
    )
    df2 = pd.concat([pd.DataFrame([pd.NA]), df2, pd.DataFrame([pd.NA])])[
        ["chrom2", "start", "end", "strand", "animal"]
    ]
    df2 = df2.astype(
        {
            "start": pd.Int64Dtype(),
            "end": pd.Int64Dtype(),
        }
    )

    assert (2, 5) == np.shape(bioframe.setdiff(df1, df1, cols1=cols1, cols2=cols1))
    assert (2, 5) == np.shape(bioframe.setdiff(df1, df2, cols1=cols1, cols2=cols2))
    assert (4, 5) == np.shape(
        bioframe.setdiff(df1, df2, on=["strand"], cols1=cols1, cols2=cols2)
    )


def test_count_overlaps():
    df1 = pd.DataFrame(
        [
            ["chr1", 8, 12, "+", "cat"],
            ["chr1", 8, 12, "-", "cat"],
            ["chrX", 1, 8, "+", "cat"],
        ],
        columns=["chrom1", "start", "end", "strand", "animal"],
    )

    df2 = pd.DataFrame(
        [
            ["chr1", 6, 10, "+", "dog"],
            ["chr1", 6, 10, "+", "dog"],
            ["chrX", 7, 10, "+", "dog"],
            ["chrX", 7, 10, "+", "dog"],
        ],
        columns=["chrom2", "start2", "end2", "strand", "animal"],
    )

    assert (
        bioframe.count_overlaps(
            df1,
            df2,
            on=None,
            cols1=("chrom1", "start", "end"),
            cols2=("chrom2", "start2", "end2"),
        )["count"].values
        == np.array([2, 2, 2])
    ).all()

    assert (
        bioframe.count_overlaps(
            df1,
            df2,
            on=["strand"],
            cols1=("chrom1", "start", "end"),
            cols2=("chrom2", "start2", "end2"),
        )["count"].values
        == np.array([2, 0, 2])
    ).all()

    assert (
        bioframe.count_overlaps(
            df1,
            df2,
            on=["strand", "animal"],
            cols1=("chrom1", "start", "end"),
            cols2=("chrom2", "start2", "end2"),
        )["count"].values
        == np.array([0, 0, 0])
    ).all()

    # overlaps with pd.NA
    counts_no_nans = bioframe.count_overlaps(
        df1,
        df2,
        on=None,
        cols1=("chrom1", "start", "end"),
        cols2=("chrom2", "start2", "end2"),
    )

    df1_na = (pd.concat([pd.DataFrame([pd.NA]), df1, pd.DataFrame([pd.NA])])).astype(
        {
            "start": pd.Int64Dtype(),
            "end": pd.Int64Dtype(),
        }
    )[["chrom1", "start", "end", "strand", "animal"]]
    df2_na = (pd.concat([pd.DataFrame([pd.NA]), df2, pd.DataFrame([pd.NA])])).astype(
        {
            "start2": pd.Int64Dtype(),
            "end2": pd.Int64Dtype(),
        }
    )[["chrom2", "start2", "end2", "strand", "animal"]]

    counts_nans_inserted_after = (
        pd.concat([pd.DataFrame([pd.NA]), counts_no_nans, pd.DataFrame([pd.NA])])
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype(),
    })[
        ["chrom1", "start", "end", "strand", "animal", "count"]
    ]

    counts_nans = bioframe.count_overlaps(
        df1_na,
        df2_na,
        on=None,
        cols1=("chrom1", "start", "end"),
        cols2=("chrom2", "start2", "end2"),
    )

    pd.testing.assert_frame_equal(
        counts_nans,
        bioframe.count_overlaps(
            df1_na,
            df2,
            on=None,
            cols1=("chrom1", "start", "end"),
            cols2=("chrom2", "start2", "end2"),
        ),
    )

    assert (counts_nans['count'].values == counts_nans_inserted_after['count'].fillna(0).values).all()



def test_pair_by_distance():
    df = pd.DataFrame(
        [
            ["chr1", 1, 3, "+", "cat"],
            ["chr1", 6, 8, "+", "skunk"],
            ["chr1", 9, 11, "-", "dog"],
        ],
        columns=["chrom", "start", "end", "strand", "animal"],
    )

    # Distance between midpoints, from_ends=False
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=4,
            min_interjacent=None,
            max_interjacent=None,
            from_ends=False,
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[6, 8, 9, 11]])
    ).all()

    # Distance between regions ends, from_ends=True
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=4,
            min_interjacent=None,
            max_interjacent=None,
            from_ends=True,
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[1, 3, 6, 8]])
    ).all()

    # Distance between midpoints is large
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=6,
            min_interjacent=None,
            max_interjacent=None,
            from_ends=False,
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[1, 3, 6, 8], [6, 8, 9, 11]])
    ).all()

    # Distance between midpoints is large
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=9,
            min_interjacent=None,
            max_interjacent=None,
            from_ends=False,
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[1, 3, 6, 8], [1, 3, 9, 11], [6, 8, 9, 11]])
    ).all()

    # Do not allow interjacent regions
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=9,
            min_interjacent=None,
            max_interjacent=0,
            from_ends=False,
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[1, 3, 6, 8], [6, 8, 9, 11]])
    ).all()

    # Strictly one interjacent region
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=9,
            min_interjacent=1,
            max_interjacent=None,
            from_ends=False,
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[1, 3, 9, 11]])
    ).all()


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
    pd.testing.assert_frame_equal(df_assigned, bioframe.assign_view(df, view_df))

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
        bioframe.assign_view(
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
        bioframe.assign_view(
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

    pd.testing.assert_frame_equal(df_sorted, bioframe.sort_bedframe(df))

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
        df_sorted, bioframe.sort_bedframe(df, view_df, view_name_col="fruit")
    )

    ### chr2 is not in the view, so if infer_assignment=False this should raise a ValueError
    assert pytest.raises(
        ValueError,
        bioframe.sort_bedframe,
        df,
        view_df,
        view_name_col="fruit",
        infer_assignment=False,
    )
