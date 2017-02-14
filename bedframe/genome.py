from __future__ import division, print_function, unicode_literals

import numpy as np
import pandas

from .util import read_table

CHROMINFO_URLS = {
    'hg38': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz',
    'hg19': 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz',
    'mm10': 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromInfo.txt.gz',
    'mm9': 'http://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/chromInfo.txt.gz',
}


def fetch_chromsizes(db, **kwargs):
    """
    Download chromosome sizes from UCSC as a ``pandas.Series``, indexed by
    chromosome label.

    """
    return read_chrominfo(CHROMINFO_URLS[db], name_index=True, **kwargs)['length']


def read_chrominfo(filepath_or_fp,
                    name_patterns=(r'^chr[0-9]+$', r'^chr[XY]$', r'^chrM$'),
                    name_index=False,
                    all_names=False,
                    **kwargs):
    """
    Parse a ``<db>.chrom.sizes`` or ``<db>.chromInfo.txt`` file from the UCSC
    database, where ``db`` is a genome assembly name.

    Input
    -----
    filepath_or_fp : str or file-like
        Path or url to text file, or buffer.
    name_patterns : sequence, optional
        Sequence of regular expressions to capture desired sequence names.
        Each corresponding set of records will be sorted in natural order.
    name_index : bool, optional
        Index table by chromosome name.
    all_names : bool, optional
        Whether to return all scaffolds listed in the file. Default is
        ``False``.

    Returns
    -------
    Data frame indexed by sequence name, with columns 'name' and 'length'.

    """
    chrom_table = read_table(filepath_or_fp, usecols=[0, 1],
                             names=['name', 'length'], **kwargs)
    if not all_names:
        parts = []
        for pattern in name_patterns:
            part = chrom_table[chrom_table['name'].str.contains(pattern)]
            part = part.iloc[argnatsort(part['name'])]
            parts.append(part)
        chrom_table = pandas.concat(parts, axis=0)
    chrom_table.insert(0, 'id', np.arange(len(chrom_table)))
    if name_index:
        chrom_table.index = chrom_table['name'].values
    return chrom_table


def binnify(chromsizes, binsize):
    """
    Divide a genome into evenly sized bins.

    Parameters
    ----------
    chromsizes : Series
        pandas Series indexed by chromosome name with chromosome lengths in bp.
    binsize : int
        size of bins in bp

    Returns
    -------
    Data frame with columns: 'chrom', 'start', 'end'.

    """
    def _each(chrom):
        clen = chromsizes.at[chrom]
        n_bins = int(np.ceil(clen / binsize))
        binedges = np.arange(0, (n_bins+1)) * binsize
        binedges[-1] = clen
        return pandas.DataFrame({
                'chrom': [chrom]*n_bins,
                'start': binedges[:-1],
                'end': binedges[1:],
            }, columns=['chrom', 'start', 'end'])
    return pandas.concat(map(_each, chromsizes.index),
                         axis=0, ignore_index=True)

