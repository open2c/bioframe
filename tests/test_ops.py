import pandas as pd
import bioframe
import pyranges as pr
import numpy as np

def bioframe_to_pyranges(df):
    pydf = df.copy()
    pydf.rename({'chrom':'Chromosome','start':'Start','end':'End'},
                  axis='columns', inplace=True)
    return pr.PyRanges(pydf)

def pyranges_overlap_to_bioframe(pydf):
    ## convert the df output by pyranges join into a bioframe-compatible format
    df = pydf.df.copy()
    df.rename({'Chromosome':'chrom_1','Start':'start_1','End':'end_1',
               'Start_b':'start_2','End_b':'end_2'},
                  axis='columns', inplace=True)
    df['chrom_1'] = df['chrom_1'].values.astype('object') #to remove categories
    df['chrom_2'] = df['chrom_1'].values    
    return df

chroms = ['chr12','chrX']
def mock_bioframe(num_entries =100):
    pos = np.random.randint(1, 1e7, size=(num_entries,2))
    df = pd.DataFrame()
    df['chrom'] = np.random.choice(chroms,num_entries)
    df['start'] = np.min(pos,axis=1)
    df['end'] = np.max(pos,axis=1)
    df.sort_values(['chrom','start'],inplace=True)
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
    assert (df1.equals(df2)==False)
    p1 = bioframe_to_pyranges(df1)
    p2 = bioframe_to_pyranges(df2)
    pp = pyranges_overlap_to_bioframe(p1.join(p2))[['chrom_1','start_1','end_1','chrom_2','start_2','end_2']]
    bb = bioframe.overlap(df1,df2)[['chrom_1','start_1','end_1','chrom_2','start_2','end_2']]
    pd.testing.assert_frame_equal(bb,pp,check_dtype=False, check_exact=True)
    print('overlap elements agree')




