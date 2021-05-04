import pandas as pd
import numpy as np
import collections
from bioframe.ops import _get_default_colnames, _verify_columns
import bioframe

def is_contained(df, view_df, view_col='parent', cols=None):
    """
    TODO: update after 
    tests if all genomic intervals in a bioframe `df` are contained in the view `view_df`.
    
    df : pandas.DataFrame
    
    view : pandas.DataFrame 
    
    <<TODO: discuss how columns from views should be handeled in general.
    << e.g. do we need a _get_default_view_colnames()
    << that gets a quadruplet(chrom,start,end,name)? >>
    
    view_col: 
    cols:
    
    Returns
    -------
    is_contained:bool
    
    """

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    df_trim = bioframe.ops.trim(df, limits=view_df, limits_region_col='name')
    is_start_trimmed = np.any(df['start'].values != df_trim['start'].values)
    is_end_trimmed   = np.any(df['end'].values   != df_trim['end'].values)

    if is_start_trimmed or is_end_trimmed:
        return False
    else:
        return True

def is_overlapping(df, cols=None):
    """
    
    tests if any genomic intervals in a bioframe `df` overlap
    
    Returns
    -------
    is_overlapping:bool
    
    """

    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    if (len(df) > len(bioframe.ops.merge(df, cols=(ck1, sk1, ek1) ))):
        return True
    else:
        return False

### discuss: needed?
# -- is_gapped(df, cols=None):

# -- is_covering(df, regions)
###     <<TODO: discuss desired behavior of complement ###

# -- is_tiling(df, regions)
# combo of: is_covering and is_contained and not is_overlapping

# -- is_sorted(df, regions=None) 
# (e.g. only sorts by chrom,start,end if no regions specified)

