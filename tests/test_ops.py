import pandas as pd
import bioframe
import pyranges as pr
import numpy as np
from io import StringIO


def bioframe_to_pyranges(df):
    pydf = df.copy()
    pydf.rename(
        {"chrom": "Chromosome", "start": "Start", "end": "End"},
        axis="columns",
        inplace=True,
    )
    return pr.PyRanges(pydf)


def pyranges_to_bioframe(pydf):
    df = pydf.df
    df.rename(
        {"Chromosome": "chrom", "Start": "start", "End": "end", "Count": "n_intervals"},
        axis="columns",
        inplace=True,
    )
    return df


def pyranges_overlap_to_bioframe(pydf):
    ## convert the df output by pyranges join into a bioframe-compatible format
    df = pydf.df.copy()
    df.rename(
        {
            "Chromosome": "chrom_1",
            "Start": "start_1",
            "End": "end_1",
            "Start_b": "start_2",
            "End_b": "end_2",
        },
        axis="columns",
        inplace=True,
    )
    df["chrom_1"] = df["chrom_1"].values.astype("object")  # to remove categories
    df["chrom_2"] = df["chrom_1"].values
    return df


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


def test_expand():
    fake_bioframe = pd.DataFrame(
        {"chrom": ["chr1", "chr1", "chr2"], "start": [1, 50, 100], "end": [5, 55, 200]}
    )
    fake_chromsizes = {"chr1": 60, "chr2": 300}
    expand_bp = 10
    fake_expanded = bioframe.expand(fake_bioframe.copy(), expand_bp, fake_chromsizes)
    print(fake_expanded)
    assert fake_expanded.iloc[0].start == 0  # don't expand below zero
    assert (
        fake_expanded.iloc[1].end == fake_chromsizes["chr1"]
    )  # don't expand above chromsize
    assert (
        fake_expanded.iloc[2].end == fake_bioframe.iloc[2].end + expand_bp
    )  # expand end normally
    assert (
        fake_expanded.iloc[2].start == fake_bioframe.iloc[2].start - expand_bp
    )  # expand start normally


def test_overlap():
    ### note does not test overlap_start or overlap_end columns of bioframe.overlap
    df1 = mock_bioframe()
    df2 = mock_bioframe()
    assert df1.equals(df2) == False
    p1 = bioframe_to_pyranges(df1)
    p2 = bioframe_to_pyranges(df2)
    pp = pyranges_overlap_to_bioframe(p1.join(p2))[
        ["chrom_1", "start_1", "end_1", "chrom_2", "start_2", "end_2"]
    ]
    bb = bioframe.overlap(df1, df2)[
        ["chrom_1", "start_1", "end_1", "chrom_2", "start_2", "end_2"]
    ]
    pd.testing.assert_frame_equal(bb, pp, check_dtype=False, check_exact=True)
    print("overlap elements agree")


def test_cluster():
    df1 = pd.DataFrame(
        [["chr1", 1, 5], ["chr1", 3, 8], ["chr1", 8, 10], ["chr1", 12, 14],],
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
    assert (
        (bioframe_to_pyranges(df1).cluster(count=True).df["Cluster"].values - 1)
        == bioframe.cluster(df1.sort_values(["chrom", "start"]))["cluster"].values
    ).all()


def test_merge():
    df1 = pd.DataFrame(
        [["chr1", 1, 5], ["chr1", 3, 8], ["chr1", 8, 10], ["chr1", 12, 14],],
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

    # test consistency with pyranges
    pd.testing.assert_frame_equal(
        pyranges_to_bioframe(bioframe_to_pyranges(df1).merge(count=True)),
        bioframe.merge(df1),
        check_dtype=False,
        check_exact=True,
    )


def test_complement():
    df1 = pd.DataFrame(
        [["chr1", 1, 5], ["chr1", 3, 8], ["chr1", 8, 10], ["chr1", 12, 14]],
        columns=["chrom", "start", "end"],
    )
    df1_chromsizes = {"chr1": 100, "chrX": 100}

    df1_complement = pd.DataFrame(
        [["chr1", 0, 1], ["chr1", 10, 12], ["chr1", 14, 100]],
        columns=["chrom", "start", "end"],
    )

    pd.testing.assert_frame_equal(
        bioframe.complement(df1, chromsizes=df1_chromsizes), df1_complement
    )

    ### test complement with two chromosomes ###
    df1.iloc[0, 0] = "chrX"
    df1_complement = pd.DataFrame(
        [
            ["chr1", 0, 3],
            ["chr1", 10, 12],
            ["chr1", 14, 100],
            ["chrX", 0, 1],
            ["chrX", 5, 100],
        ],
        columns=["chrom", "start", "end"],
    )
    pd.testing.assert_frame_equal(
        bioframe.complement(df1, chromsizes=df1_chromsizes), df1_complement
    )


def test_closest():
    df1 = pd.DataFrame([["chr1", 1, 5],], columns=["chrom", "start", "end"])

    df2 = pd.DataFrame(
        [["chr1", 4, 8], ["chr1", 10, 11]], columns=["chrom", "start", "end"]
    )

    ### closest(df1,df2,k=1) ###
    d = """chrom_1  start_1  end_1 chrom_2  start_2  end_2  distance
        0    chr1        1      5    chr1        4      8         0"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, bioframe.closest(df1, df2, k=1))

    ### closest(df1,df2, ignore_overlaps=True)) ###
    d = """chrom_1 start_1 end_1   chrom_2 start_2 end_2   distance
        0   chr1    1   5   chr1    10  11  5"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, bioframe.closest(df1, df2, ignore_overlaps=True))

    ### closest(df1,df2,k=2) ###
    d = """chrom_1 start_1 end_1   chrom_2 start_2 end_2   distance
            0   chr1    1   5   chr1    4   8   0
            1   chr1    1   5   chr1    10  11  5"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, bioframe.closest(df1, df2, k=2))

    ### closest(df2,df1) ###
    d = """chrom_1  start_1 end_1   chrom_2 start_2 end_2   distance
            0   chr1    4   8   chr1    1   5   0
            1   chr1    10  11  chr1    1   5   5 """
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, bioframe.closest(df2, df1))

    ### change first interval to new chrom ###
    df2.iloc[0, 0] = "chrA"
    d = """chrom_1 start_1 end_1   chrom_2 start_2 end_2   distance
              0   chr1    1   5   chr1    10  11  5"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, bioframe.closest(df1, df2, k=1))


def test_coverage():


    #### coverage does not exceed length of original interval
    df1 = pd.DataFrame([
        ['chr1', 3, 8]],
        columns=['chrom', 'start', 'end']
    )
    df2 = pd.DataFrame([
        ['chr1', 2, 10] ],
        columns=['chrom', 'start', 'end']
    )
    d = """chrom    start   end coverage    n_overlaps
         0  chr1    3   8   5   1"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, bioframe.coverage(df1, df2))


    ### coverage of interval on different chrom returns zero for coverage and n_overlaps 
    df1 = pd.DataFrame([
        ['chr1', 3, 8]],
        columns=['chrom', 'start', 'end']
    )
    df2 = pd.DataFrame([
        ['chrX', 3, 8] ],
        columns=['chrom', 'start', 'end']
    )
    d = """chrom    start   end coverage    n_overlaps
         0  chr1    3   8   0   0"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, bioframe.coverage(df1, df2))


    ### currently does not pass, as n_overlaps not calculated correctly 
    ### when a second overlap starts within the first
    df1 = pd.DataFrame([["chr1", 3, 8]], columns=["chrom", "start", "end"])
    df2 = pd.DataFrame(
        [["chr1", 3, 5], ["chr1", 5, 8]], columns=["chrom", "start", "end"]
    )

    d = """chrom    start   end coverage    n_overlaps
         0  chr1    3   8   5   2"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(df, bioframe.coverage(df1, df2))



