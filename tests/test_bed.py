import os
import tempfile

import pandas as pd
import pytest

import bioframe


def test_involution():
    with tempfile.TemporaryDirectory() as directory:
        for schema in ['narrowPeak', 'bed12']:
            bf = bioframe.read_table(f'tests/test_data/{schema}.bed',
                                     schema=schema)
            fname = os.path.join(directory, f'{schema}.bed')
            bioframe.to_bed(bf, fname)
            involution = bioframe.read_table(fname, schema=schema)
            pd.testing.assert_frame_equal(bf, involution)


def test_chrom_validators():
    with tempfile.TemporaryDirectory() as directory:
        bf = bioframe.read_table('tests/test_data/bed12.bed', schema='bed12')
        bf.loc[0, 'chrom'] = 'value with space'
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        bf.loc[0, 'chrom'] = '' # must be non empty
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        bf.loc[0, 'chrom'] = 'a'*300 # must be shorter than 256
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))


def test_end_validators():
    with tempfile.TemporaryDirectory() as directory:
        bf = bioframe.read_table('tests/test_data/bed12.bed', schema='bed12')
        bf.loc[0, 'end'] = 10 # end must be after start
        bf.loc[0, 'start'] = 11
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))


def test_name_validators():
    with tempfile.TemporaryDirectory() as directory:
        bf = bioframe.read_table('tests/test_data/bed12.bed', schema='bed12')
        bf.loc[0, 'name'] = '' # must not be empty
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        bf.loc[0, 'name'] = 'a'*300 # must be less than 255 char
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))


def test_score_validators():
    with tempfile.TemporaryDirectory() as directory:
        bf = bioframe.read_table('tests/test_data/bed12.bed', schema='bed12')
        # negative value is enforced by the normal types

        bf.loc[0, 'score'] = 1001
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'), strict_score=True)

        bf.loc[0, 'score'] = '.' # enforced to be a number by the types
        with pytest.raises(TypeError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))


def test_strand_validators():
    with tempfile.TemporaryDirectory() as directory:
        bf = bioframe.read_table('tests/test_data/bed12.bed', schema='bed12')
        bf.loc[0, 'strand'] = '*'
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))


def test_thick_validators():
    with tempfile.TemporaryDirectory() as directory:
        for direction in ['Start', 'End']:
            bf = bioframe.read_table('tests/test_data/bed12.bed', schema='bed12')
            bf.loc[0, 'start'] = 100
            bf.loc[0, 'end'] = 1000
            bf.loc[0, f'thick{direction}'] = 1001
            with pytest.raises(ValueError):
                bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

            bf.loc[0, f'thick{direction}'] = 99
            with pytest.raises(ValueError):
                bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))


def test_itemRgb_validators():
    with tempfile.TemporaryDirectory() as directory:
        bf = bioframe.read_table('tests/test_data/bed12.bed', schema='bed12')
        bf["itemRgb"] = bf["itemRgb"].astype(str)
        bf.loc[0, 'itemRgb'] = 'a,12,13' # must be integers
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        bf.loc[0, 'itemRgb'] = '12,13' # must be 1 or 3 integers
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        bf.loc[0, 'itemRgb'] = '12,13,14,15' # must be 1 or 3 integers
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        bf.loc[0, 'itemRgb'] = '12,13,300' # must be between 0 and 255
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        bf.loc[0, 'itemRgb'] = '300' # must be between 0 and 255
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))


def test_blockCount_validators():
    with tempfile.TemporaryDirectory() as directory:
        bf = bioframe.read_table('tests/test_data/bed12.bed', schema='bed12')
        bf.loc[0, 'blockCount'] = 0
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))


def test_blockSizes_validators():
    with tempfile.TemporaryDirectory() as directory:
        bf = bioframe.read_table('tests/test_data/bed12.bed', schema='bed12')
        bf.loc[0, 'blockCount'] = 2
        bf.loc[0, 'blockSizes'] = '2,a,'
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        bf.loc[0, 'blockCount'] = 2
        bf.loc[0, 'blockSizes'] = '2,2,2,'
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))


def test_blockStarts_validators():
    with tempfile.TemporaryDirectory() as directory:
        bf = bioframe.read_table('tests/test_data/bed12.bed', schema='bed12')
        bf.loc[0, 'blockCount'] = 2
        bf.loc[0, 'blockSizes'] = '2,4,'
        bf.loc[0, 'blockStarts'] = '0,a,'
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        bf.loc[0, 'blockCount'] = 2
        bf.loc[0, 'blockSizes'] = '1,1,'
        bf.loc[0, 'blockStarts'] = '0,2,5,'
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        # ends after end
        bf.loc[0, 'start'] = 1
        bf.loc[0, 'end'] = 10
        bf.loc[0, 'blockCount'] = 1
        bf.loc[0, 'blockSizes'] = '100,'
        bf.loc[0, 'blockStarts'] = '0,'
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        # ends before end
        bf.loc[0, 'start'] = 1
        bf.loc[0, 'end'] = 10
        bf.loc[0, 'blockCount'] = 1
        bf.loc[0, 'blockSizes'] = '1,'
        bf.loc[0, 'blockStarts'] = '0,'
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))

        # overlap
        bf.loc[0, 'start'] = 1
        bf.loc[0, 'end'] = 10
        bf.loc[0, 'blockCount'] = 2
        bf.loc[0, 'blockSizes'] = '5,5,'
        bf.loc[0, 'blockStarts'] = '0,1,'
        with pytest.raises(ValueError):
            bioframe.to_bed(bf, os.path.join(directory, 'foo.bed'))
