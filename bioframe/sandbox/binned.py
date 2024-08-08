import numpy as np
import pandas as pd

import scipy.signal


from ..core import arrops


def stackup_track(track_df, windows, window_size, value_col, fill_val=np.nan):
    """
    Extract equisized windows from a track dataframe and stack them up.

    Parameters
    ----------
    track_df : pandas.DataFrame
        DataFrame with a column for the values of the track and a column for the genome position.
    windows : pandas.DataFrame
        A set of genomic intervals used to extract the windows. Note that only the midpoints of the intervals are used!
    window_size : int
        Size of the windows.
    value_col : str
        Name of the column in track_df that contains the values of the track.
    fill_val : int or float, optional
        Value to use to fill the end of a partially-filled window. Defaults to np.nan.

    Returns
    -------
    stack : np.array
        A 2D array of shape (n_windows, window_size) containing the stacked track windows.
    """

    bin_widths = (track_df['end'] - track_df['start'])[
        (track_df['chrom'] == track_df['chrom'].shift(-1))
        ].values
    resolution = bin_widths[0]
    if not np.isclose(bin_widths, resolution).all():
        raise ValueError('Bin widths are not uniform')

    window_size_pixels = window_size // resolution
    n_windows = len(windows)

    out_stack = np.full((n_windows, window_size_pixels), fill_val, dtype=np.float_)

    track_gb = track_df.groupby('chrom')
    regions_gb = windows.reset_index(drop=True).groupby('chrom') # there could be more elegant solutions not involving reset_index, i.e. windows.index[group]
    groups = set(track_gb.groups.keys()) & set(regions_gb.groups.keys())

    for chrom in groups:
        track_df_chrom = track_gb.get_group(chrom)
        windows_df_chrom = regions_gb.get_group(chrom)

        track_bin_mids_chrom = (track_df_chrom.start.values + resolution // 2)
        window_mids = (windows_df_chrom.start.values + windows_df_chrom.end.values) // 2

        mids_idxs = np.searchsorted(track_bin_mids_chrom, window_mids, side='left')
        start_idxs = mids_idxs - window_size_pixels // 2
        end_idxs = start_idxs + window_size_pixels 

        track_chrom = track_df_chrom[value_col].values
        padded_track_chrom = np.pad(track_chrom, window_size_pixels, mode='constant', constant_values=fill_val)

        stack_chrom = padded_track_chrom[arrops.arange_multi(start_idxs, end_idxs) + window_size_pixels # shift by the padding width
                                         ].reshape(-1, window_size_pixels)

        out_stack[regions_gb.indices[chrom]] = stack_chrom

    return out_stack


def smooth_track(df, value_col, window, min_periods, weight_col='weight', suffix='_smoothed', win_type=None, sum_kwargs={}):
    """
    Smooths a track in a DataFrame using a rolling window.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the track data.
        value_col (str): The name of the column to be smoothed.
        window (int): The size of the rolling window.
        min_periods (int): The minimum number of non-null values required to compute a smoothed value.
        weight_col (str, optional): The name of the column used for weighting. Defaults to 'weight'.
        suffix (str, optional): The suffix to be added to the smoothed column names. Defaults to '_smoothed'.
        win_type (str, optional): The window type to be used for smoothing. Defaults to None.
        sum_kwargs (dict, optional): Additional keyword arguments to be passed to the sum function. Defaults to {}.

    Returns:
        pandas.DataFrame: The DataFrame with the smoothed track.
    """


    gb=df.groupby('chrom')
    dfs = []
    for chrom in df.chrom.unique():
        dfc = pd.DataFrame(gb.get_group(chrom))
        dfc[weight_col] = 1 - dfc[value_col].isnull().astype(np.float_)
        for col in (value_col, weight_col):
            dfc[col+suffix] = dfc[col].rolling(
                window=window, 
                min_periods=min_periods,
                win_type=win_type,
                center=True).sum(**sum_kwargs)
        dfc[value_col+suffix] /= dfc[weight_col + suffix]
        dfs.append(dfc)

    return pd.concat(dfs)


def smooth_gauss(df, value_col, sigma, min_frac=0.2, window_sigma=6):
    """
    Smooths a DataFrame column using a Gaussian window.

    Args:
        df (pandas.DataFrame): The DataFrame containing the data.
        value_col (str): The name of the column to be smoothed.
        sigma (float): The standard deviation of the Gaussian window.
        min_frac (float, optional): The minimum fraction of non-null values required to compute a smoothed value. Defaults to 0.2.
        window_sigma (float, optional): The size of the window in terms of standard deviations. Defaults to 6.

    Returns:
        pandas.DataFrame: The DataFrame with the smoothed column.

    """
    df = smooth_track(
        df, 
        value_col, 
        window_sigma*sigma, 
        int(np.round(window_sigma*sigma*min_frac)), 
        suffix='_smooth_gaussian',
        win_type='gaussian', 
        sum_kwargs={'std':sigma})
    return df


def smooth_square(df, value_col, window, min_frac=0.2):
    """
    Smooths the values in a DataFrame column using a square window.

    Args:
        df (pandas.DataFrame): The DataFrame containing the data.
        value_col (str): The name of the column to be smoothed.
        window (int): The size of the square window for smoothing.
        min_frac (float, optional): The minimum fraction of non-null values required to compute a smoothed value. Defaults to 0.2.

    Returns:
        pandas.DataFrame: The DataFrame with the smoothed values in a new column.

    """
    df = smooth_track(
        df, 
        value_col, 
        window, 
        int(np.round(window*min_frac)), 
        suffix='_smooth_square')
    return df


def find_peaks(track_df, value_col):
    """
    Find peaks in a genomic track and calculate their prominence.

    Parameters:
    track_df (pandas.DataFrame): a bedgraph-like genomic dataframe.
    value_col (str): the column in 'track_df' that contains the values for which the peak prominences will be calculated.

    Returns:
    track_df (pandas.DataFrame): The modified DataFrame with a new column named 'value_col'+'prominence', which contains the calculated peak prominences.
    """

    track_gb = track_df.groupby('chrom')

    proms = np.full(len(track_df), np.nan, dtype=np.float_)

    for chrom in track_gb.groups:
        track_df_chrom = track_gb.get_group(chrom)
        value_track_chrom = track_df_chrom[value_col].values

        peaks_chrom = scipy.signal.find_peaks(value_track_chrom)[0]
        proms_chrom = scipy.signal.peak_prominences(value_track_chrom, peaks_chrom)[0]
        
        proms_track_chrom = np.full(len(track_df_chrom), np.nan, dtype=np.float_)
        proms_track_chrom[peaks_chrom] = proms_chrom
        proms[track_gb.indices[chrom]] = proms_track_chrom

    track_df[value_col+'_prominence'] = proms

    return track_df


def select_bins(track_df, features, value_col):
    """
    Selects the values from a track DataFrame corresponding to the midpoints of features.
    It uses a more efficient algorithm than overlap, but it requires the features to be sorted.

    Args:
        track_df (pandas.DataFrame): DataFrame containing the track data.
        features (pandas.DataFrame): DataFrame containing the features data.
        value_col (str): Name of the column in track_df containing the values.

    Returns:
        numpy.ndarray: Array of selected values from the track data.
    """
    
    n_features = len(features)

    out_signal = np.full((n_features,), np.nan, dtype=np.float_)

    track_gb = track_df.groupby('chrom')
    features_gb = features.reset_index(drop=True).groupby('chrom')
    groups = set(track_gb.groups.keys()) & set(features_gb.groups.keys())

    for chrom in groups:
        track_df_chrom = track_gb.get_group(chrom)
        features_df_chrom = features_gb.get_group(chrom)

        features_mids = (features_df_chrom.start.values + features_df_chrom.end.values) // 2

        mids_idxs = np.searchsorted(track_df_chrom.start.values, features_mids, side='left')
        # add a check that the mids are within the bins

        out_signal_chrom = track_df_chrom[value_col].values[mids_idxs]

        out_signal[features_gb.indices[chrom]] = out_signal_chrom

    return out_signal