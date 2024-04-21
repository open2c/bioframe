import os.path as op

import numpy as np
import pandas as pd
import pytest

import bioframe

testdir = op.realpath(op.dirname(__file__))


def test_make_chromarms():
    ### test the case where columns have different names
    df = pd.DataFrame(
        [["chrX", 0, 8]],
        columns=["chromosome", "lo", "hi"],
    )
    mids = pd.DataFrame([["chrX", 4]], columns=["chromosome", "loc"])
    arms = pd.DataFrame(
        [
            ["chrX", 0, 4, "chrX_p"],
            ["chrX", 4, 8, "chrX_q"],
        ],
        columns=["chrom", "start", "end", "name"],
    )
    arms = arms.astype({"start": pd.Int64Dtype(), "end": pd.Int64Dtype()})

    # test passing 3 columns
    result = bioframe.make_chromarms(
        df,
        mids,
        cols_chroms=["chromosome", "lo", "hi"],
        cols_mids=["chromosome", "loc"],
    )
    pd.testing.assert_frame_equal(
        result, arms.rename(columns={"chrom": "chromosome", "start": "lo", "end": "hi"})
    )

    # test passing 2 columns
    result = bioframe.make_chromarms(
        df,
        mids,
        cols_chroms=["chromosome", "hi"],
        cols_mids=["chromosome", "loc"],
    )
    pd.testing.assert_frame_equal(
        result,
        arms.rename(columns={"chrom": "chromosome"}),
    )

    # test for passing Series or dict
    result = bioframe.make_chromarms(
        pd.Series({"chrX": 8}), mids, cols_mids=["chromosome", "loc"]
    )
    pd.testing.assert_frame_equal(arms, result)

    result = bioframe.make_chromarms(pd.Series({"chrX": 8}), pd.Series({"chrX": 4}))
    pd.testing.assert_frame_equal(arms, result)

    bioframe.make_chromarms({"chrX": 8}, mids, cols_mids=["chromosome", "loc"])
    pd.testing.assert_frame_equal(arms, result)

    bioframe.make_chromarms({"chrX": 8}, pd.Series({"chrX": 4}))
    pd.testing.assert_frame_equal(arms, result)

    bioframe.make_chromarms({"chrX": 8}, {"chrX": 4})
    pd.testing.assert_frame_equal(arms, result)


def test_binnify():
    chromsizes = bioframe.read_chromsizes(
        testdir + "/test_data/test.chrom.sizes", filter_chroms=False
    )
    assert len(chromsizes) == 2
    assert len(bioframe.binnify(chromsizes, int(np.max(chromsizes.values)))) == len(
        chromsizes
    )
    assert len(bioframe.binnify(chromsizes, int(np.min(chromsizes.values)))) == (
        len(chromsizes) + 1
    )
    assert len(bioframe.binnify(chromsizes, 1)) == np.sum(chromsizes.values)


def test_digest():
    pytest.importorskip("Bio")
    fasta_records = bioframe.load_fasta(testdir + "/test_data/test.fa")
    assert len(fasta_records) == 2
    ### no HindIII sites in the test.fa fasta records, so shouldn't change shape[0]
    assert bioframe.digest(fasta_records, "HindIII").shape == (2, 3)
    ### one DpnII site on chrTEST2, shape[0] should increase by one
    assert bioframe.digest(fasta_records, "DpnII").shape == (3, 3)
    ### DpnII site is on chrTEST2 position 3, first interval of chrTEST2 should end at 3
    assert bioframe.digest(fasta_records, "DpnII").iloc[1].end == 3


def test_frac_mapped():
    pytest.importorskip("pysam")
    chromsizes = bioframe.read_chromsizes(
        testdir + "/test_data/test.chrom.sizes", filter_chroms=False
    )
    fasta_records = bioframe.load_fasta(testdir + "/test_data/test.fa")

    unmapped = np.array([1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
    assert (
        unmapped
        == bioframe.frac_mapped(
            bioframe.binnify(chromsizes, 1), fasta_records, return_input=False
        ).values
    ).all()

    unmapped = np.array([0.8, 0.8, 0])
    assert (
        unmapped
        == bioframe.frac_mapped(
            bioframe.binnify(chromsizes, 5), fasta_records, return_input=False
        ).values
    ).all()

    unmapped = np.array([0.8, 4 / 7])
    assert (
        unmapped
        == bioframe.frac_mapped(
            bioframe.binnify(chromsizes, 7), fasta_records, return_input=False
        ).values
    ).all()


def test_frac_gc():
    pytest.importorskip("pysam")
    chromsizes = bioframe.read_chromsizes(
        testdir + "/test_data/test.chrom.sizes", filter_chroms=False
    )
    fasta_records = bioframe.load_fasta(testdir + "/test_data/test.fa")

    unmapped_bp = (
        0
        == bioframe.frac_mapped(
            bioframe.binnify(chromsizes, 1), fasta_records, return_input=False
        ).values
    )
    assert np.isnan(
        bioframe.frac_gc(
            bioframe.binnify(chromsizes, 1),
            fasta_records,
            return_input=False,
            mapped_only=True,
        ).values[unmapped_bp]
    ).all()

    ## mapped_only=True should ignore N or return np.nan if interval only contains N
    np.testing.assert_equal(
        np.array([0.5, 0.5, np.nan]),
        bioframe.frac_gc(
            bioframe.binnify(chromsizes, 5),
            fasta_records,
            return_input=False,
            mapped_only=True,
        ).values,
    )

    assert (
        np.array([0.5, 0.5])
        == bioframe.frac_gc(
            bioframe.binnify(chromsizes, 7),
            fasta_records,
            return_input=False,
            mapped_only=True,
        ).values
    ).all()

    ## mapped_only=False should count N as zero
    assert (
        np.array([0.4, 0.4, 0])
        == bioframe.frac_gc(
            bioframe.binnify(chromsizes, 5),
            fasta_records,
            return_input=False,
            mapped_only=False,
        ).values
    ).all()

    assert (
        np.array([0.4, 2 / 7])
        == bioframe.frac_gc(
            bioframe.binnify(chromsizes, 7),
            fasta_records,
            return_input=False,
            mapped_only=False,
        ).values
    ).all()


def test_seq_gc():
    assert 0 == bioframe.seq_gc("AT")
    assert np.isnan(bioframe.seq_gc("NNN"))
    assert 1 == bioframe.seq_gc("NGnC")
    assert 0.5 == bioframe.seq_gc("GTCA")
    assert 0.25 == bioframe.seq_gc("nnnNgTCa", mapped_only=False)
    with pytest.raises(ValueError):
        bioframe.seq_gc(["A", "T"])
    with pytest.raises(ValueError):
        bioframe.seq_gc(np.array("ATGC"))


### todo: test frac_gene_coverage(bintable, mrna):
### currently broken


def test_pair_by_distance():
    df = pd.DataFrame(
        [
            ["chr1", 1, 3, "+", "cat"],
            ["chr1", 6, 8, "+", "skunk"],
            ["chr1", 9, 11, "-", "dog"],
        ],
        columns=["chrom", "start", "end", "strand", "animal"],
    )

    # Distance between midpoints
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=4,
            min_intervening=None,
            max_intervening=None,
            relative_to="midpoints",
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[6, 8, 9, 11]])
    ).all()

    # Distance between regions endpoints
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=4,
            min_intervening=None,
            max_intervening=None,
            relative_to="endpoints",
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[1, 3, 6, 8]])
    ).all()

    # Distance between midpoints is large
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=6,
            min_intervening=None,
            max_intervening=None,
            relative_to="midpoints",
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[1, 3, 6, 8], [6, 8, 9, 11]])
    ).all()

    # Distance between midpoints is large
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=9,
            min_intervening=None,
            max_intervening=None,
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[1, 3, 6, 8], [1, 3, 9, 11], [6, 8, 9, 11]])
    ).all()

    # Do not allow intervening regions
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=9,
            min_intervening=None,
            max_intervening=0,
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[1, 3, 6, 8], [6, 8, 9, 11]])
    ).all()

    # Strictly one intervening region
    assert (
        bioframe.pair_by_distance(
            df,
            min_sep=1,
            max_sep=9,
            min_intervening=1,
            max_intervening=None,
        )[["start_1", "end_1", "start_2", "end_2"]].values
        == np.array([[1, 3, 9, 11]])
    ).all()

    # no negative min_sep
    with pytest.raises(ValueError):
        bioframe.pair_by_distance(df, min_sep=-1, max_sep=9)

    # no min_sep > max_sep
    with pytest.raises(ValueError):
        bioframe.pair_by_distance(df, min_sep=12, max_sep=9)

    # no min_intervening > max_intervening
    with pytest.raises(ValueError):
        bioframe.pair_by_distance(
            df, min_sep=0, max_sep=9, min_intervening=10, max_intervening=9
        )
