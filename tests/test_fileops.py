from io import StringIO

import pandas as pd
import pytest

import bioframe


############# tests #####################
def test_read_table():
    d = """chr1\nchr2\nchr2"""
    assert bioframe.read_table(StringIO(d), schema="bed3").shape == (3, 3)

    # raise a value error if any columns are filled with all NA
    with pytest.raises(ValueError):
        bioframe.read_table(StringIO(d), schema="bed3", schema_is_strict=True)

    # fill with nans to appropriate size if schema_is_strict=False (aka the default)
    d = """chr1      5    10
           chr1     10   20
           chr2    30  40"""
    assert bioframe.read_table(StringIO(d), schema="bed3", sep="\\s+").shape == (3, 3)
    assert bioframe.read_table(StringIO(d), schema="bed6", sep="\\s+").shape == (3, 6)
    assert bioframe.read_table(StringIO(d), schema="bed12", sep="\\s+").shape == (3, 12)

    # bedpe has 10 columns
    d = """chr1    5    10  chr2   5   10   interval1  .  +  -
           chr1    10   20  chr1   5   10   interval2  .  +  -
           chr2    30   40  chr2   5   10   interval3  12  +  -
        """
    assert bioframe.read_table(
        StringIO(d), schema="bedpe", sep=r"\s+", schema_is_strict=True
    ).shape == (3, 10)


def test_read_chromsizes():
    d = """chr1\nchr2\nchr2"""
    with pytest.raises(ValueError):
        bioframe.read_chromsizes(StringIO(d))

    d = """chr1\t1\nchr3\t2\nchr2\t3\n """
    chromsizes = bioframe.read_chromsizes(StringIO(d))
    assert isinstance(chromsizes, pd.Series)
    assert chromsizes.name == "length"
    assert list(chromsizes.index) == ["chr1", "chr2", "chr3"]
    assert list(chromsizes.values) == [1, 3, 2]


def test_read_beds():
    # Checking that we properly read common bed schemas
    schemas = ['narrowPeak', 'jaspar', 'bed9', 'bed12']

    for schema in schemas:
        _ = bioframe.read_table(f'tests/test_data/{schema}.bed', schema=schema,
                                schema_is_strict=True)


def test_read_sam():
    # SAM file taken from https://github.com/samtools/samtools/blob/develop/examples/toy.sam
    _ = bioframe.read_alignment('tests/test_data/toy.sam')


def test_read_bam():
    # converted toy.sam via `samtools view -bS toy.sam > toy.bam;
    # index file created with `samtools index toy.bam`
    _ = bioframe.read_alignment('tests/test_data/toy.bam')
