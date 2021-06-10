from io import StringIO

import pandas as pd
import numpy as np
import pytest

from bioframe.core.construction import from_any
from bioframe.core import construction


def test_add_ucsc_name_column():
    df = pd.DataFrame(
        {"chrom": [f"chr{i}" for i in range(3)], "start": [1, 2, 3], "end": [4, 5, 6]}
    )

    pd.testing.assert_series_equal(
        construction.add_ucsc_name_column(df)["name"],
        pd.Series(
            data=["chr0:1-4", "chr1:2-5", "chr2:3-6"], index=[0, 1, 2], name="name"
        ),
    )


def test_any():

    ### tests copied from old parse_regions
    # main functionality: convert to dataframe and create name
    df = pd.DataFrame(
        {"chrom": [f"chr{i}" for i in range(3)], "start": [1, 2, 3], "end": [4, 5, 6]}
    )
    parsed = from_any(df, names_as_ucsc=True, name_col="regions")
    assert parsed.iloc[0]["chrom"] == "chr0"
    assert parsed.iloc[0]["regions"] == "chr0:1-4"

    # re-create dataframe from UCSC name alone
    df2 = pd.DataFrame({"regions": parsed["regions"].values})
    assert (from_any(df2, names_as_ucsc=True, name_col="regions") == parsed).all().all()

    # re-parsing results yields the same
    assert (from_any(parsed) == parsed).all().all()

    # None or False will be parsed
    assert from_any([("chr1", None, 5)], fill_null={"chr1": 10})["start"].values[0] == 0

    # pull end from chromsizes
    p2 = from_any([("chr1", 5, None)], fill_null={"chr1": 40}, names_as_ucsc=True)
    assert list(p2.values[0]) == ["chr1", 5, 40, "chr1:5-40"]

    # We could keep things as None if chromsizes were not proviced
    p3 = from_any(["chr1", "chr2"], fill_null=False)
    assert list(p3.values[0]) == ["chr1", None, None, "chr1"]

    # we can force UCSC names
    p4 = from_any([("chr1", 0, 5)], names_as_ucsc=True)
    assert list(p4.values[0]) == ["chr1", 0, 5, "chr1:0-5"]

    # nothing happens: name was autocompleted
    p5 = from_any([("chr1", 0, 5)], names_as_ucsc=False)
    assert list(p5.values[0]) == ["chr1", 0, 5, "chr1"]

    # forcing UCSC names
    p7 = from_any([("chr1")], fill_null={"chr1": 5}, names_as_ucsc=True)
    assert list(p7.values[0]) == ["chr1", 0, 5, "chr1:0-5"]

    # kept the strange name
    p8 = from_any(["chr1:1,000,000-4M"])
    assert list(p8.values[0]) == ["chr1", 1000000, 4000000, "chr1:1,000,000-4M"]

    p9 = from_any(["chr1"])
    assert list(p9.values[0]) == ["chr1", None, None, "chr1"]

    with pytest.raises(ValueError):
        from_any([("ch1", 1, 2, "chr1:1-2", "puppies")])  # puppies are not allowed

    with pytest.raises(ValueError):
        from_any([("chr1", 5, None)], fill_null={"chr2": 40})


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
        sanitized_df1, construction.sanitize_bedframe(df1, drop_null=True)
    )

    # keep rows with null, but recast
    sanitized_df1 = df1.astype(
        {"chrom": str, "start": pd.Int64Dtype(), "end": pd.Int64Dtype()}
    )
    pd.testing.assert_frame_equal(sanitized_df1, construction.sanitize_bedframe(df1))

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
        construction.sanitize_bedframe(
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
    assert construction.sanitize_bedframe(
        df1, start_exceed_end_action="drop", drop_null=True
    ).empty


def test_make_viewframe():

    # test dict input
    d = """            chrom  start  end        name
    0    chrTESTX      0   10    chrTESTX
    1  chrTESTX_p      0   12  chrTESTX_p"""
    view_df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(
        view_df.copy(),
        construction.make_viewframe({"chrTESTX": 10, "chrTESTX_p": 12}),
    )

    # test list input
    region_list = [("chrTESTX", 0, 10), ("chrTESTX_p", 0, 12)]
    pd.testing.assert_frame_equal(
        view_df.copy(),
        construction.make_viewframe(region_list),
    )

    # test pd.Series input
    chromsizes = pd.Series(data=[5, 8], index=["chrTESTXq", "chrTEST_2p"])
    d = """      chrom  start  end        name
    0  chrTESTXq      0    5   chrTESTXq
    1   chrTEST_2p      0    8  chrTEST_2p"""
    view_df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(
        view_df.copy(), construction.make_viewframe(chromsizes)
    )

    d = """          chrom   start   end name
    0   chrTESTXq   0   5   chrTESTXq:0-5
    1   chrTEST_2p  0   8   chrTEST_2p:0-8"""
    view_df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(
        view_df.copy(),
        construction.make_viewframe(chromsizes, view_names_as_ucsc=True),
    )

    # test pd.DataFrame input
    pd.testing.assert_frame_equal(view_df.copy(), construction.make_viewframe(view_df))
