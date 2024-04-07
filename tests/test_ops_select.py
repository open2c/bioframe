import pandas as pd
import pytest

import bioframe


def test_select():
    df = pd.DataFrame(
        [["chrX", 3, 8], ["chr1", 4, 5], ["chrX", 1, 5]],
        columns=["chrom", "start", "end"],
    )

    result = pd.DataFrame([["chr1", 4, 5]], columns=["chrom", "start", "end"])
    pd.testing.assert_frame_equal(
        result, bioframe.select(df, "chr1:4-10").reset_index(drop=True)
    )

    result = pd.DataFrame(
        [["chrX", 3, 8], ["chrX", 1, 5]], columns=["chrom", "start", "end"]
    )
    pd.testing.assert_frame_equal(
        result, bioframe.select(df, "chrX").reset_index(drop=True)
    )

    result = pd.DataFrame(
        [["chrX", 3, 8], ["chrX", 1, 5]], columns=["chrom", "start", "end"]
    )
    pd.testing.assert_frame_equal(
        result, bioframe.select(df, "chrX:4-6").reset_index(drop=True)
    )

    # Query range not in the dataframe
    assert len(bioframe.select(df, "chrZ")) == 0
    assert len(bioframe.select(df, "chr1:100-1000")) == 0
    assert len(bioframe.select(df, "chr1:1-3")) == 0

    # Invalid query range
    with pytest.raises(ValueError):
        bioframe.select(df, "chr1:1-0")


def test_select__with_colnames():
    ### select with non-standard column names
    new_names = ["chr", "chrstart", "chrend"]
    df = pd.DataFrame(
        [["chrX", 3, 8], ["chr1", 4, 5], ["chrX", 1, 5]],
        columns=new_names,
    )
    result = pd.DataFrame(
        [["chrX", 3, 8], ["chrX", 1, 5]],
        columns=new_names,
    )
    pd.testing.assert_frame_equal(
        result, bioframe.select(df, "chrX:4-6", cols=new_names).reset_index(drop=True)
    )
    pd.testing.assert_frame_equal(
        result, bioframe.select(df, "chrX", cols=new_names).reset_index(drop=True)
    )


def test_select__with_nulls():
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

    result = pd.DataFrame(
        [["chr1", -6, 12, "chr1p"]],
        columns=colnames,
    ).astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    pd.testing.assert_frame_equal(
        result, bioframe.select(df, "chr1:0-1").reset_index(drop=True)
    )


def test_select__mask_indices_labels():
    df = pd.DataFrame(
        [["chrX", 3, 8], ["chr1", 4, 5], ["chrX", 1, 5]],
        columns=["chrom", "start", "end"],
    )

    region = "chr1:4-10"
    answer = pd.DataFrame([["chr1", 4, 5]], columns=["chrom", "start", "end"])

    result = bioframe.select(df, region)
    pd.testing.assert_frame_equal(answer, result.reset_index(drop=True))
    mask = bioframe.select_mask(df, region)
    pd.testing.assert_frame_equal(answer, df.loc[mask].reset_index(drop=True))
    labels = bioframe.select_labels(df, region)
    pd.testing.assert_frame_equal(answer, df.loc[labels].reset_index(drop=True))
    idx = bioframe.select_indices(df, region)
    pd.testing.assert_frame_equal(answer, df.iloc[idx].reset_index(drop=True))


def test_select__query_intervals_are_half_open():
    df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2"],
            "start": [0, 10, 10, 20, 30, 40, 50, 60],
            "end": [10, 20, 20, 30, 40, 50, 60, 70],
            "name": ["a", "b", "A", "B", "C", "D", "E", "F"],
        }
    )

    result = bioframe.select(df, "chr1")
    assert (result["name"] == ["a", "b"]).all()

    result = bioframe.select(df, "chr2:20-70")
    assert (result["name"] == ["B", "C", "D", "E", "F"]).all()

    result = bioframe.select(df, "chr2:20-75")
    assert (result["name"] == ["B", "C", "D", "E", "F"]).all()

    result = bioframe.select(df, "chr2:20-")
    assert (result.index == [3, 4, 5, 6, 7]).all()

    result = bioframe.select(df, "chr2:20-30")
    assert (result["name"] == ["B"]).all()

    result = bioframe.select(df, "chr2:20-40")
    assert (result["name"] == ["B", "C"]).all()

    result = bioframe.select(df, "chr2:20-45")
    assert (result["name"] == ["B", "C", "D"]).all()

    result = bioframe.select(df, "chr2:19-45")
    assert (result["name"] == ["A", "B", "C", "D"]).all()

    result = bioframe.select(df, "chr2:25-45")
    assert (result["name"] == ["B", "C", "D"]).all()

    result = bioframe.select(df, "chr2:25-50")
    assert (result["name"] == ["B", "C", "D"]).all()

    result = bioframe.select(df, "chr2:25-51")
    assert (result["name"] == ["B", "C", "D", "E"]).all()


def test_select__with_point_intervals():
    # Dataframe containing "point intervals"
    df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2"],
            "start": [0, 10, 10, 20, 30, 40, 50, 60],
            "end": [10, 10, 20, 30, 40, 50, 50, 70],
            "name": ["a", "b", "A", "B", "C", "D", "E", "F"],
        }
    )
    result = bioframe.select(df, "chr1")
    assert (result["name"] == ["a", "b"]).all()

    result = bioframe.select(df, "chr1:4-10")
    assert (result["name"] == ["a"]).all()

    result = bioframe.select(df, "chr1:4-4")
    assert (result["name"] == ["a"]).all()

    result = bioframe.select(df, "chr1:10-15")
    assert (result["name"] == ["b"]).all()

    result = bioframe.select(df, "chr2:20-70")
    assert (result["name"] == ["B", "C", "D", "E", "F"]).all()

    result = bioframe.select(df, "chr2:49-70")
    assert (result["name"] == ["D", "E", "F"]).all()

    result = bioframe.select(df, "chr2:50-70")
    assert (result["name"] == ["E", "F"]).all()

    result = bioframe.select(df, "chr2:50-51")
    assert (result["name"] == ["E"]).all()

    result = bioframe.select(df, "chr2:50-50")
    assert (result["name"] == ["E"]).all()


def test_select__with_points():
    # Dataframe of points
    df = pd.DataFrame(
        [["chrX", 3, "A"], ["chr1", 4, "C"], ["chrX", 1, "B"]],
        columns=["chrom", "pos", "name"],
    )

    result = bioframe.select(df, "chr1:4-10", cols=["chrom", "pos", "pos"])
    assert (result["name"] == ["C"]).all()

    result = bioframe.select(df, "chr1:3-10", cols=["chrom", "pos", "pos"])
    assert (result["name"] == ["C"]).all()

    result = bioframe.select(df, "chr1:4-4", cols=["chrom", "pos", "pos"])
    assert (result["name"] == ["C"]).all()
