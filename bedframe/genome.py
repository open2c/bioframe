from __future__ import division, print_function, unicode_literals

import numpy as np
import pandas

from .util import read_table
from .io import read_chrominfo

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
        clen = chromsizes[chrom]
        n_bins = int(np.ceil(clen / binsize))
        binidxs = np.arange(0, (n_bins+1))
        binedges = binidxs * binsize
        binedges[-1] = clen
        return pandas.DataFrame({
                'chrom': [chrom]*n_bins,
                'start': binedges[:-1],
                'end': binedges[1:],
                'binidx' : binidxs[:-1],
            }, columns=['chrom', 'start', 'end','binidx',])

    chromTable = pandas.concat(map(_each, chromsizes.keys()),
                               axis=0, ignore_index=True)

#    chromTable['chrom'] = pandas.Categorical(
#        chromTable.chrom, 
#        categories=list(chromsizes.keys()), 
#        ordered=True)
    return chromTable


def digest(fasta_records, enzyme):
    # http://biopython.org/DIST/docs/cookbook/Restriction.html#mozTocId447698
    chroms = fasta_records.keys()
    try:
        cut_finder = getattr(biorst, enzyme).search
    except AttributeError:
        raise ValueError('Unknown enzyme name: {}'.format(enzyme))

    def _each(chrom):
        seq = bioseq.Seq(str(fasta_records[chrom]))
        cuts = np.r_[0, np.array(cut_finder(seq)) + 1, len(seq)].astype(int)
        n_frags = len(cuts) - 1
        frags = pandas.DataFrame({
            'chrom': [chrom] * n_frags,
            'start': cuts[:-1],
            'end': cuts[1:],
            'binidx': np.arange(0, n_frags)},
            columns=['chrom', 'start', 'end', 'binidx'])
        return frags

    chromTable = pandas.concat(map(_each, chroms), axis=0, ignore_index=True)
    chromTable['chrom'] = pandas.Categorical(chromTable.chrom, categories=chroms, ordered=True)
    return chromTable


def frac_mapped(bintable, fasta_records):
    def _each(bin):
        s = str(fasta_records[bin.chrom][bin.start:bin.end])
        nbases = len(s)
        n = s.count('N')
        n += s.count('n')
        return (nbases - n) / nbases
    return bintable.apply(_each, axis=1)


def frac_gc(bintable, fasta_records, mapped_only=True):
    def _each(bin):
        s = str(fasta_records[bin.chrom][bin.start:bin.end])
        g = s.count('G')
        g += s.count('g')
        c = s.count('C')
        c += s.count('c')
        nbases = len(s)
        if mapped_only:
            n = s.count('N')
            n += s.count('n')
            nbases -= n
        return (g + c) / nbases if nbases > 0 else np.nan
    return bintable.apply(_each, axis=1)

class Genome:
    '''
    Tasks: 
    
    '''
    def __init__(self,
            chromsizes='',
            centromeres=None,
            bins=None,
            fasta_path=None,
            mapped_frac=False,
            GC=False,
            ):
        pass
    
    def from_cache(path):
        pass
    
    @property
    def chroms(self)
        return self._chroms
    
    @property
    def chromarms(self)
        return self._chromarms
 
    @property
    def bins(self)
        return self._bins
