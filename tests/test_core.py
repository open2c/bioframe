
import pandas as pd
import numpy as np
from io import StringIO
import bioframe.core
import bioframe
import os.path as op

testdir = op.realpath(op.dirname(__file__))

#def test_is_valid():

#def test_sanitize():

#def test_is_view():

def test_make_view():
    # test dict input

    # test pd.Series input
    chromsizes = pd.Series(data=[5,8], index=['chrTESTXq','chrTEST_2p'])
    d = """      chrom  start  end        name
    0  chrTESTX      0    5   chrTESTXq
    1   chrTEST      0    8  chrTEST_2p"""
    view_df = pd.read_csv(StringIO(d), sep=r"\s+")
    pd.testing.assert_frame_equal(view_df.copy(), bioframe.core.make_view(chromsizes))

    # test pd.DataFrame input
    pd.testing.assert_frame_equal(view_df.copy(), bioframe.core.make_view(view_df))


