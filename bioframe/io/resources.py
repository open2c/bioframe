from __future__ import division, print_function
from .formats import read_chromsizes, read_gapfile


def fetch_chromsizes(db, **kwargs):
    """
    Download chromosome sizes from UCSC as a ``pandas.Series``, indexed by
    chromosome label.

    """
    return read_chromsizes(
        'http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/chromInfo.txt.gz'.format(db), 
        **kwargs)


def fetch_gaps(db, **kwargs):
    return read_gapfile(
    	'http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/gap.txt.gz'.format(db),
    	**kwargs)
