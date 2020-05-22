import collections

import numpy as np
import pandas as pd

from . import arrops
from . import ops
from ._region import parse_region

import six


def make_chromarms(chromsizes, mids, binsize=None, suffixes=("p", "q")):
    """
    Split chromosomes into chromosome arms

    Parameters
    ----------
    chromsizes : pandas.Series
        Series mapping chromosomes to lengths in bp.

    mids : dict-like
        Mapping of chromosomes to midpoint locations.

    binsize : int, optional
        Round midpoints to nearest bin edge for compatibility with a given
        bin grid.

    suffixes : tuple, optional
        Suffixes to name chromosome arms. Defaults to p and q.

    Returns
    -------
    4-column BED-like DataFrame (chrom, start, end, name).
    Arm names are chromosome names + suffix.
    Any chromosome not included in ``mids`` will be omitted.

    """
    chromosomes = [chrom for chrom in chromsizes.index if chrom in mids]

    p_arms = [[chrom, 0, mids[chrom], chrom + suffixes[0]] for chrom in chromosomes]
    if binsize is not None:
        for x in p_arms:
            x[2] = int(round(x[2] / binsize)) * binsize

    q_arms = [
        [chrom, mids[chrom], chromsizes[chrom], chrom + suffixes[1]]
        for chrom in chromosomes
    ]
    if binsize is not None:
        for x in q_arms:
            x[1] = int(round(x[1] / binsize)) * binsize

    interleaved = [*sum(zip(p_arms, q_arms), ())]

    return pd.DataFrame(interleaved, columns=["chrom", "start", "end", "name"])


def binnify(chromsizes, binsize, rel_ids=False):
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
    bintable : pandas.DataFrame with columns: 'chrom', 'start', 'end'.

    """

    if type(binsize) is not int:
        raise ValueError("binsize must be int")

    def _each(chrom):
        clen = chromsizes[chrom]
        n_bins = int(np.ceil(clen / binsize))
        binedges = np.arange(0, (n_bins + 1)) * binsize
        binedges[-1] = clen
        return pd.DataFrame(
            {"chrom": [chrom] * n_bins, "start": binedges[:-1], "end": binedges[1:]},
            columns=["chrom", "start", "end"],
        )

    bintable = pd.concat(map(_each, chromsizes.keys()), axis=0, ignore_index=True)

    if rel_ids:
        bintable["rel_id"] = bintable.groupby("chrom").cumcount()

    # if as_cat:
    #     bintable['chrom'] = pd.Categorical(
    #         bintable['chrom'],
    #         categories=list(chromsizes.keys()),
    #         ordered=True)

    return bintable


def digest(fasta_records, enzyme):
    """
    Divide a genome into restriction fragments.

    Parameters
    ----------
    fasta_records : OrderedDict
        Dictionary of chromosome names to sequence records.
        Created by: bioframe.load_fasta('/path/to/fasta.fa')

    enzyme: str
        Name of restriction enzyme.

    Returns
    -------
    Dataframe with columns: 'chrom', 'start', 'end'.

    """
    import Bio.Restriction as biorst
    import Bio.Seq as bioseq

    # http://biopython.org/DIST/docs/cookbook/Restriction.html#mozTocId447698
    chroms = fasta_records.keys()
    try:
        cut_finder = getattr(biorst, enzyme).search
    except AttributeError:
        raise ValueError("Unknown enzyme name: {}".format(enzyme))

    def _each(chrom):
        seq = bioseq.Seq(str(fasta_records[chrom][:]))
        cuts = np.r_[0, np.array(cut_finder(seq)) + 1, len(seq)].astype(int)
        n_frags = len(cuts) - 1

        frags = pd.DataFrame(
            {"chrom": [chrom] * n_frags, "start": cuts[:-1], "end": cuts[1:]},
            columns=["chrom", "start", "end"],
        )
        return frags

    return pd.concat(map(_each, chroms), axis=0, ignore_index=True)


def frac_mapped(df, fasta_records, return_input=True):
    """
    Calculate the fraction of mapped base-pairs for each interval in a dataframe. 

    Parameters
    ----------
    df : pandas.DataFrame
        A sets of genomic intervals stored as a DataFrame.

    fasta_records : OrderedDict
        Dictionary of chromosome names to sequence records.
        Created by: bioframe.load_fasta('/path/to/fasta.fa')

    return_input: bool
        if False, only return Series named frac_mapped.

    Returns
    -------
    df_mapped : pd.DataFrame 
        Original dataframe with new column 'frac_mapped' appended.

    """

    if not set(df["chrom"].values).issubset(set(fasta_records.keys())):
        return ValueError(
            "chrom from intervals not in fasta_records: double-check genome agreement"
        )

    def _each(bin):
        s = str(fasta_records[bin.chrom][bin.start : bin.end])
        nbases = len(s)
        n = s.count("N")
        n += s.count("n")
        return (nbases - n) / nbases

    if return_input:
        return pd.concat(
            [df, df.apply(_each, axis=1).rename("frac_mapped", inplace=True)],
            axis="columns",
        )
    else:
        return df.apply(_each, axis=1).rename("frac_mapped", inplace=True)


def frac_gc(df, fasta_records, mapped_only=True, return_input=True):
    """
    Calculate the fraction of GC basepairs for each interval in a dataframe. 

    Parameters
    ----------
    df : pandas.DataFrame
        A sets of genomic intervals stored as a DataFrame.

    fasta_records : OrderedDict
        Dictionary of chromosome names to sequence records.
        Created by: bioframe.load_fasta('/path/to/fasta.fa')

    mapped_only: bool
        if True, ignore 'N' in the fasta_records for calculation. 
        if True and there are no mapped base-pairs in an interval, return np.nan.

    return_input: bool
        if False, only return Series named frac_mapped.

    Returns
    -------
    df_mapped : pd.DataFrame 
        Original dataframe with new column 'frac_mapped' appended.
    
    """
    if not set(df["chrom"].values).issubset(set(fasta_records.keys())):
        return ValueError(
            "chrom from intervals not in fasta_records: double-check genome agreement"
        )

    def _each(chrom_group):
        chrom = chrom_group.name
        seq = fasta_records[chrom]
        gc = []
        for _, bin in chrom_group.iterrows():
            s = str(seq[bin.start : bin.end])
            g = s.count("G")
            g += s.count("g")
            c = s.count("C")
            c += s.count("c")
            nbases = len(s)
            if mapped_only:
                n = s.count("N")
                n += s.count("n")
                nbases -= n
            gc.append((g + c) / nbases if nbases > 0 else np.nan)
        return gc

    out = df.groupby("chrom", sort=False).apply(_each)

    if return_input:
        return pd.concat(
            [df, pd.Series(data=np.concatenate(out), index=df.index).rename("GC")],
            axis="columns",
        )
    else:
        return pd.Series(data=np.concatenate(out), index=df.index).rename("GC")


def frac_gene_coverage(bintable, mrna):
    """
    Calculate number and fraction of overlaps by genes for a set of intervals.

    Parameters
    ----------
    bintable : set of intervals
    mrna: str
        Name of genome.

    Returns
    -------
    XXXX

    """

    raise ValueError("implementation currently broken!")

    if isinstance(mrna, six.string_types):
        from .resources import UCSCClient

        mrna = (
            UCSCClient(mrna)
            .fetch_mrna()
            .rename(columns={"tName": "chrom", "tStart": "start", "tEnd": "end"})
        )

    #### currently broken ####
    bintable = ops.coverage(
        bintable,
        mrna,
        out={"input": "input", "count": "gene_count", "coverage": "gene_coverage"},
    )

    return bintable
