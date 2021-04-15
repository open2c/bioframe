import pandas as pd
import numpy as np
import collections
from bioframe.ops import _get_default_colnames, _verify_columns

def _verify_column_dtypes(df, cols = None):
    """
    Checks that a dataframe has chrom, start, end columns and with valid dtypes.
    
    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.

    """
    
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    _verify_columns(df, [ck1, sk1, ek1])

    column_dtypes = df.dtypes
    if not  np.any(
            [pd.api.types.is_string_dtype(column_dtypes[ck1]),
             pd.api.types.is_object_dtype(column_dtypes[ck1]),
             pd.api.types.is_categorical_dtype(column_dtypes[ck1])]):
        raise TypeError("invalid df['chrom'] dtype, must be object, string, or categorical")
    if not pd.api.types.is_integer_dtype(column_dtypes[sk1]):
        raise TypeError("invalid df['start'] dtype, must be integer")
    if not pd.api.types.is_integer_dtype(column_dtypes[ek1]):
        raise TypeError("invalid df['end'] dtype, must be integer")

def _verify_genomic_interval_df(df, drop_invalid = False, flip_invalid = False, cols = None): 
    """
    Checks that df has chrom, start, end columns, that their dtypes are valid,
    and that all starts < ends. 

    drop_invalid:bool
        returns a copy of the dataframe without invalid invervals. default False.
        
    flip_invalid:bool
        flips intervals where start<end and returns a copy of the dataframe where these intervals have been flipped.
        default False.
        
    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.
        
    """
    
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    _verify_column_dtypes(df)

    out_df = df.copy()
    if drop_invalid:
        out_df = out_df.iloc[((df[ek1] - df[sk1]) >= 0).values]
    elif ((df[ek1] - df[sk1]) <0).any(): 
        if flip_invalid:
            inds = ((df[ek1] - df[sk1]) < 0).values
            out_df.loc[inds,[sk1,ek1]] = out_df.loc[inds,[ek1,sk1]].values
        else:
            raise ValueError("Invalid genomic interval dataframe: starts exceed ends for "+str(np.sum(((df[ek1] - df[sk1]) <0)))+" intervals")
    return out_df

def _verify_region_df(region_df, region_name_col='name', cols=None):
    """
    Checks that region_df has chrom, start, end, region_name_col
    and that entries in the region_name_col are unique.
    """
    
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols
    
    _verify_columns(region_df, [ck1, sk1, ek1, region_name_col])
    
    _verify_genomic_interval_df(region_df)

    if pd.isna(region_df).values.any():
        raise ValueError('Invalid region dataframe: cannot contain NAs')
    
    if (len(set(region_df[region_name_col])) < len(region_df[region_name_col].values)) :
        raise ValueError('Invalid region dataframe: entries in region_df[region_name_col] must be unique')

    return region_df
            

def make_region_df(regions, infer_chroms_from_regions=True, region_name_col='name', cols=None):
    """
    Makes and validates a region_df, where supported input types for regions are:
    - a dictionary where keys are strings and values are integers {str:int},
    specifying regions (chrom, 0, end, name)
    - a dictionary where keys are strings and values are tuples of integers {str:(int,int)},
    specifying regions (chrom, start, end, name)
    - a dataFrame, skips to validation step
    
    infer_chroms_from_regions:bool
        attemps to strip 'p' or 'q' from chrom string. if False, region_name_col specifies chrom as well.
        default True.
        
    region_name_col:str
        specifies 
        
    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.
        
    
    """
    ck1, sk1, ek1 = _get_default_colnames() if cols is None else cols

    if type(regions) is dict:
        data = []
        for k, v in dict(regions).items():
            name = k
            if infer_chroms_from_regions:
                chrom = k.split('_')[0].replace('p','').replace('q','')
            else: 
                chrom = k
            if isinstance(v, (tuple, list, np.ndarray)):
                start = v[0]
                end = v[1]
            elif np.isscalar(v):
                start = 0
                end = v
            else:
                raise ValueError("Unknown dict format: {type(v)}")
            data.append([chrom, start, end, name])

        regions_df = pd.DataFrame(data, columns=[ck1, sk1, ek1, region_name_col])
        
    elif type(regions) is pd.core.frame.DataFrame:
        regions_df = regions.copy()
        
    else: 
        raise ValueError("Unknown region type: {type(v)}")  
        
    _verify_region_df(regions_df, region_name_col= region_name_col, cols=(ck1,sk1,ek1))
    
    return regions_df
