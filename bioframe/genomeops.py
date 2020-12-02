import numpy as np
import pandas as pd

from . import ops


def make_chromarms(
    chromsizes,
    midpoints,
    cols_chroms=("chrom", "length"),
    cols_mids=("chrom", "mids"),
    suffixes=("_p", "_q"),
):
    """
    Split chromosomes into chromosome arms

    Parameters
    ----------
    chromsizes : pandas.Dataframe or pandas.Series
        If pandas.Series, a map from chromosomes to lengths in bp.
        If pandas.Dataframe, a dataframe with columns 'chrom' and 'length'.

    midpoints : pandas.Dataframe or dict-like
        Mapping of chromosomes to midpoint locations.

    suffixes : tuple, optional
        Suffixes to name chromosome arms. Defaults to p and q.

    Returns
    -------
    df_chromarms
        4-column BED-like DataFrame (chrom, start, end, name).
        Arm names are chromosome names + suffix.
        Any chromosome not included in ``mids`` will be omitted.

    """
    ck1, sk1 = cols_chroms
    ck2, sk2 = cols_mids

    if isinstance(chromsizes, pd.Series):
        df_chroms = (
            pd.DataFrame(chromsizes).reset_index().rename(columns={"index": ck1})
        )
    elif isinstance(chromsizes, pd.DataFrame):
        df_chroms = chromsizes.copy()
    else:
        raise ValueError("unknown input type for chromsizes")

    if isinstance(midpoints, dict):
        df_mids = pd.DataFrame.from_dict(midpoints, orient="index", columns=[sk2])
        df_mids.reset_index(inplace=True)
        df_mids.rename(columns={"index": ck2}, inplace=True)
    elif isinstance(midpoints, pd.DataFrame):
        df_mids = midpoints.copy()

    ops._verify_columns(df_mids, [ck2, sk2])
    ops._verify_columns(df_chroms, [ck1, sk1])

    df_chroms["start"] = 0
    df_chroms["end"] = df_chroms[sk1].values

    df_chromarms = ops.split(
        df_chroms,
        df_mids,
        add_names=True,
        cols=(ck1, "start", "end"),
        cols_points=(ck2, sk2),
        suffixes=suffixes,
    )
    df_chromarms["name"].replace(r"[\:\[].*?[\)\_]", "", regex=True, inplace=True)
    df_chromarms.drop(columns=["index_2", "length"], inplace=True)
    return df_chromarms


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
    try:
        import Bio.Restriction as biorst
        import Bio.Seq as bioseq
    except ImportError:
        raise ImportError("Biopython is required to use digest")

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
        return (nbases - n) / nbases if nbases > 0 else 0

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


def frac_gene_coverage(df, ucsc_mrna):
    """
    Calculate number and fraction of overlaps by genes for a set of intervals stored in a dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Set of genomic intervals stored as a dataframe.

    ucsc_mrna: str or DataFrame
        Name of UCSC genome or all_mrna.txt dataframe from UCSC or similar.

    Returns
    -------
    df_gene_coverage : pd.DataFrame

    """
    if isinstance(ucsc_mrna, str):
        from .io.resources import UCSCClient

        mrna = UCSCClient(ucsc_mrna).fetch_mrna()
    else:
        mrna = ucsc_mrna

    mrna = mrna.rename(columns={"tName": "chrom", "tStart": "start", "tEnd": "end"})
    df_gene_coverage = ops.coverage(df, mrna)
    df_gene_coverage = ops.count_overlaps(df_gene_coverage, mrna)

    return df_gene_coverage
