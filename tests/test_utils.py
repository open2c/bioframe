from io import StringIO
import os.path as op

import pandas as pd
import numpy as np
import pytest

import bioframe

testdir = op.realpath(op.dirname(__file__))

### todo: test make_chromarms(chromsizes, mids, binsize=None, suffixes=("p", "q")):


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


### todo: test frac_gene_coverage(bintable, mrna):
### currently broken
