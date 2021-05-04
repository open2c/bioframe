import pandas as pd
import bioframe
import pyranges as pr
import numpy as np
from io import StringIO
import bioframe.definitions


#todo: def test_is_contained():


def test_is_overlapping():


    d = """chrom  start  end
         0  chr1      3    6
         1  chr1     5   10
         2  chr2    5  10"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")
    #print(df)

    assert( bioframe.definitions.is_overlapping(df) is True)
    
    d = """chrom  start  end
         0  chr1      3    6
         1  chr1     7   10
         2  chr2    5  10"""
    df = pd.read_csv(StringIO(d), sep=r"\s+")

    assert( bioframe.definitions.is_overlapping(df) is False)

