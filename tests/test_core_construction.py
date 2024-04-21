from io import StringIO

import numpy as np
import pandas as pd
import pytest

from bioframe.core import construction
from bioframe.core.construction import from_any


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
    parsed = from_any(df)
    assert "name" not in parsed.columns
    assert parsed.iloc[0]["chrom"] == "chr0"

    # re-create dataframe from UCSC name alone
    df2 = pd.DataFrame(
        {
            "regions": construction.add_ucsc_name_column(parsed, name_col="regions")[
                "regions"
            ].values
        }
    )
    assert (
        (from_any(df2, name_col="regions")[["chrom", "start", "end"]] == parsed)
        .all()
        .all()
    )

    # re-parsing results yields the same
    assert (from_any(parsed) == parsed).all().all()

    # extra columns don't get overwritten
    df["name"] = "test-value"
    assert (from_any(df)["name"] == df["name"]).all()

    # None or False will be parsed
    assert from_any([("chr1", None, 5)], fill_null={"chr1": 10})["start"].values[0] == 0

    # pull end from chromsizes
    p2 = from_any([("chr1", 5, None)], fill_null={"chr1": 40})
    assert list(p2.values[0]) == ["chr1", 5, 40]

    # We could keep things as None if chromsizes were not proviced
    p3 = from_any(["chr1", "chr2"], fill_null=False)
    assert list(p3.values[0]) == ["chr1", None, None]

    # parse the strange name
    p8 = from_any(["chr1:1,000,000-4M"])
    assert list(p8.values[0]) == ["chr1", 1000000, 4000000]

    p9 = from_any(["chr1"])
    assert list(p9.values[0]) == ["chr1", None, None]

    with pytest.raises(ValueError):
        from_any([("ch1", 1, 2, "chr1:1-2", "puppies")])  # puppies are not allowed

    with pytest.raises(ValueError):
        from_any([("chr1", 5, None)], fill_null={"chr2": 40})

    # input tuple of tuples
    p2 = from_any((("chr1", 5, 10), ("chrX", 10, 20)))
    assert list(p2.values[0]) == ["chr1", 5, 10]

    # input tuple of lists
    p2 = from_any((["chr1", 5, 10], ["chrX", 10, 20]))
    assert list(p2.values[0]) == ["chr1", 5, 10]

    # input tuple of ucsc strings
    p2 = from_any(("chr1:5-10",))
    assert list(p2.values[0]) == ["chr1", 5, 10]

    # input single tuple
    p2 = from_any(("chr1", 5, 10))
    assert list(p2.values[0]) == ["chr1", 5, 10]


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
    sanitized_df1 = pd.DataFrame(
        [
            ["chr1", 10, 20],
            ["chr1", 10, 20],
            [pd.NA, pd.NA, pd.NA],
            [pd.NA, pd.NA, pd.NA],
        ],
        columns=["chrom", "start", "end"],
    )
    sanitized_df1 = sanitized_df1.astype(
        {"chrom": object, "start": pd.Int64Dtype(), "end": pd.Int64Dtype()}
    )
    pd.testing.assert_frame_equal(
        sanitized_df1.fillna(-1), construction.sanitize_bedframe(df1).fillna(-1)
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
    view_df = pd.DataFrame(
        [
            ["chrTESTX", 0, 10, "chrTESTX:0-10"],
            ["chrTESTX_p", 0, 12, "chrTESTX_p:0-12"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    pd.testing.assert_frame_equal(
        view_df.copy(),
        construction.make_viewframe(
            {"chrTESTX": 10, "chrTESTX_p": 12}, name_style="ucsc"
        ),
    )

    # test list input
    region_list = [("chrTESTX", 0, 10), ("chrTESTX_p", 0, 12)]
    pd.testing.assert_frame_equal(
        view_df.copy(),
        construction.make_viewframe(region_list, name_style="ucsc"),
    )

    # test pd.Series input
    chromsizes = pd.Series(data=[5, 8], index=["chrTESTXq", "chrTEST_2p"])
    d = """      chrom  start  end        name
    0  chrTESTXq      0    5   chrTESTXq
    1   chrTEST_2p      0    8  chrTEST_2p"""
    view_df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(
        view_df.copy(), construction.make_viewframe(chromsizes, name_style=None)
    )

    d = """          chrom   start   end name
    0   chrTESTXq   0   5   chrTESTXq:0-5
    1   chrTEST_2p  0   8   chrTEST_2p:0-8"""
    view_df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(
        view_df.copy(),
        construction.make_viewframe(chromsizes, name_style="UCSC"),
    )

    # test pd.DataFrame input
    pd.testing.assert_frame_equal(view_df.copy(), construction.make_viewframe(view_df))

    # if you provide unique names, this is accepted unchanged by make_viewframe
    view_df = pd.DataFrame(
        [["chrTESTX", 0, 10, "chrTEST_1"], ["chrTESTY", 0, 12, "chrTEST_2"]],
        columns=["chrom", "start", "end", "name"],
    )

    region_list = [("chrTESTX", 0, 10, "chrTEST_1"), ("chrTESTY", 0, 12, "chrTEST_2")]

    pd.testing.assert_frame_equal(
        view_df.copy(), construction.make_viewframe(region_list)
    )

    pd.testing.assert_frame_equal(view_df.copy(), construction.make_viewframe(view_df))

    pd.testing.assert_frame_equal(
        view_df.copy(),
        construction.make_viewframe(
            view_df, check_bounds={"chrTESTX": 11, "chrTESTY": 13}
        ),
    )

    with pytest.raises(ValueError):
        construction.make_viewframe(
            view_df, check_bounds={"chrTESTX": 9, "chrTESTY": 13}
        )
