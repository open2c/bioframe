import numpy as np
import pandas as pd
import collections

from . import ops
from .core.specs import _get_default_colnames, _verify_columns

__all__ = [
    "make_chromarms",
    "binnify",
    "digest",
    "frac_mapped",
    "frac_gc",
    "seq_gc",
    "frac_gene_coverage",
    "pair_by_distance",
]


def make_chromarms(
    chromsizes,
    midpoints,
    cols_chroms=("chrom", "length"),
    cols_mids=("chrom", "mid"),
    suffixes=("_p", "_q"),
):
    """
    Split chromosomes into chromosome arms.

    Parameters
    ----------
    chromsizes : pandas.Dataframe or pandas.Series
        If pandas.Series, a map from chromosomes to lengths in bp.
        If pandas.Dataframe, a dataframe with columns defined by cols_chroms.
        If cols_chroms is a triplet (e.g. 'chrom','start','end'), then
        values in chromsizes[cols_chroms[1]].values must all be zero.

    midpoints : pandas.Dataframe or dict-like
        Mapping of chromosomes to midpoint (aka centromere) locations.
        If pandas.Series, a map from chromosomes to midpoints in bp.
        If pandas.Dataframe, a dataframe with columns defined by cols_mids.

    cols_chroms : (str, str) or (str, str, str)
        Two columns

    suffixes : tuple, optional
        Suffixes to name chromosome arms. Defaults to p and q.

    Returns
    -------
    df_chromarms
        4-column BED-like DataFrame (chrom, start, end, name).
        Arm names are chromosome names + suffix.
        Any chromosome not included in ``mids`` will be not be split.

    """
    columns_to_drop = ["index", "sub_index_"]
    if len(cols_chroms) == 2:
        ck1, sk1 = cols_chroms
    elif len(cols_chroms) == 3:
        ck1, sk1, ek1 = cols_chroms

    if isinstance(chromsizes, pd.Series):
        df_chroms = (
            pd.DataFrame(chromsizes).reset_index().rename(columns={"index": ck1})
        )
    elif isinstance(chromsizes, pd.DataFrame):
        df_chroms = chromsizes.copy()
    else:
        raise ValueError("unknown input type for chromsizes")

    if len(cols_chroms) == 2:
        _verify_columns(df_chroms, [ck1, sk1])
        columns_to_drop += [sk1]
        df_chroms["end"] = df_chroms[sk1].values
        df_chroms["start"] = 0
        sk1, ek1 = "start", "end"
    elif len(cols_chroms) == 3:
        ck1, sk1, ek1 = cols_chroms
        _verify_columns(df_chroms, [ck1, sk1, ek1], unique_cols=True)
        if any((df_chroms[sk1].values != 0)):
            raise ValueError("all values in starts column must be zero")
    else:
        raise ValueError("invalid number of cols_chroms")

    ck2, sk2 = cols_mids
    if isinstance(midpoints, dict):
        df_mids = pd.DataFrame.from_dict(midpoints, orient="index", columns=[sk2])
        df_mids.reset_index(inplace=True)
        df_mids.rename(columns={"index": ck2}, inplace=True)
    elif isinstance(midpoints, pd.DataFrame):
        df_mids = midpoints.copy()
    else:
        raise ValueError("unknown input type for midpoints")
    _verify_columns(df_mids, [ck2, sk2])
    df_mids["start"] = df_mids[sk2]
    df_mids["end"] = df_mids[sk2]

    df_chromarms = ops.subtract(
        df_chroms,
        df_mids,
        cols1=(ck1, sk1, ek1),
        cols2=(ck2, "start", "end"),
        return_index=True,
    )
    if df_chromarms["sub_index_"].max() > 1:
        raise ValueError(
            "chromosome split into more than two arms, double-check midpoints"
        )
    df_chromarms["name"] = df_chromarms[ck1] + [
        suffixes[i] for i in df_chromarms["sub_index_"].values
    ]
    # df_chromarms.drop(columns=columns_to_drop, inplace=True)
    return df_chromarms[[ck1, sk1, ek1, "name"]]


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

    if not isinstance(binsize, int):
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
    if not isinstance(fasta_records, dict):
        raise ValueError(
            "fasta records must be provided as an OrderedDict, can be created by bioframe.load_fasta"
        )
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
        raise ValueError(
            "chrom from intervals not in fasta_records: double-check genome agreement"
        )
    if not isinstance(fasta_records, dict):
        raise ValueError(
            "fasta records must be provided as an OrderedDict, can be created by bioframe.load_fasta"
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
        Original dataframe with new column 'GC' appended.

    """
    if not set(df["chrom"].values).issubset(set(fasta_records.keys())):
        raise ValueError(
            "chrom from intervals not in fasta_records: double-check genome agreement"
        )
    if not isinstance(fasta_records, dict):
        raise ValueError(
            "fasta records must be provided as an OrderedDict, can be created by bioframe.load_fasta"
        )

    def _each(chrom_group):
        chrom = chrom_group.name
        seq = fasta_records[chrom]
        seq = str(seq[:])
        gc = []
        for _, bin in chrom_group.iterrows():
            s = seq[bin.start : bin.end]
            gc.append(seq_gc(s, mapped_only=mapped_only))
        return gc

    out = df.groupby("chrom", sort=False).apply(_each)

    if return_input:
        return pd.concat(
            [df, pd.Series(data=np.concatenate(out), index=df.index).rename("GC")],
            axis="columns",
        )
    else:
        return pd.Series(data=np.concatenate(out), index=df.index).rename("GC")


def seq_gc(seq, mapped_only=True):
    """
    Calculate the fraction of GC basepairs for a string of nucleotides.

    Parameters
    ----------
    seq : str
        Basepair input

    mapped_only: bool
        if True, ignore 'N' in the sequence for calculation.
        if True and there are no mapped base-pairs, return np.nan.

    Returns
    -------
    gc : float
        calculated gc content.

    """
    if not type(seq) == str:
        raise ValueError("reformat input sequence as a str")
    g = seq.count("G")
    g += seq.count("g")
    c = seq.count("C")
    c += seq.count("c")
    nbases = len(seq)
    if mapped_only:
        n = seq.count("N")
        n += seq.count("n")
        nbases -= n
    return (g + c) / nbases if nbases > 0 else np.nan


def frac_gene_coverage(df, ucsc_mrna):
    """
    Calculate number and fraction of overlaps by predicted and verified
    RNA isoforms for a set of intervals stored in a dataframe.

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


def pair_by_distance(
    df,
    min_sep,
    max_sep,
    min_intervening=None,
    max_intervening=None,
    relative_to="midpoints",
    cols=None,
    suffixes=("_1", "_2"),
):
    """
    From a dataframe of genomic intervals, find all unique pairs of intervals
    that are between ``min_sep`` and ``max_sep`` bp separated from each other.

    Parameters
    ----------
    df : pandas.DataFrame
        A BED-like dataframe.
    min_sep, max_sep : int
        Minimum and maximum separation between intervals in bp.
        Min > 0 and Max >= Min.
    min_intervening, max_intervening : int
        Minimum and maximum number of intervening intervals separating pairs.
        Min > 0 and Max >= Min.
    relative_to : str,
        Whether to calculate distances between interval "midpoints" or "endpoints".
        Default "midpoints".
    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals, provided separately for each set. The default
        values are 'chrom', 'start', 'end'.
    suffixes : (str, str), optional
        The column name suffixes for the two interval sets in the output.
        The first interval of each output pair is always upstream of the
        second.

    Returns
    -------
    pandas.DataFrame
        A BEDPE-like dataframe of paired intervals from ``df``.

    """
    ck, sk, ek = _get_default_colnames() if cols is None else cols

    # Sort intervals by genomic coordinates
    df = df.sort_values([ck, sk, ek]).reset_index(drop=True)

    if min_sep >= max_sep:
        raise ValueError("min_sep must be less than max_sep")
    if min_sep < 0:
        raise ValueError("min_sep must be >=0")

    if min_intervening is None:
        min_intervening = 0
    if max_intervening is None:
        max_intervening = df.index.max()
    if min_intervening > max_intervening:
        raise ValueError("min_intervening must be less or equal to max_intervening")
    if min_intervening < 0:
        raise ValueError("min_intervening must be >=0")

    mids = (df[sk] + df[ek]) // 2

    # For each interval, generate a probe interval on its right
    if relative_to == "endpoints":
        print("endpoint")
        ref = df[ek]
    elif relative_to == "midpoints":
        ref = mids
    else:
        raise ValueError("relative_to must either specify 'midpoints' or 'endpoints' ")
    right_probe = df[[ck]].copy()
    right_probe[sk] = ref + min_sep // 2
    right_probe[ek] = ref + (max_sep + 1) // 2

    # For each interval, also generate a probe interval on its left
    if relative_to == "endpoints":
        ref = df[sk]
    elif relative_to == "midpoints":
        ref = mids
    else:
        raise ValueError("relative_to must either specify 'midpoints' or 'endpoints' ")
    left_probe = df[[ck]].copy()
    left_probe[sk] = ref - max_sep // 2
    left_probe[ek] = ref - (min_sep + 1) // 2

    # Intersect right-handed probes (from intervals on the left)
    # with left-handed probes (from intervals on the right)
    idxs = ops.overlap(
        right_probe,
        left_probe,
        suffixes=suffixes,
        how="inner",
        return_index=True,
        return_input=False,
    )

    # Select only the pairs that are separated by
    # at least min_intervening intervals and no more than max_intervening intervals
    idxs["intervening"] = (
        np.abs(idxs[f"index{suffixes[0]}"] - idxs[f"index{suffixes[1]}"]) - 1
    )
    idxs = idxs.query(
        f"intervening<={max_intervening} and intervening>={min_intervening}"
    )

    left_ivals = (
        df.iloc[idxs[f"index{suffixes[0]}"].values]
        .rename(columns=lambda x: x + suffixes[0])
        .reset_index(drop=True)
    )
    right_ivals = (
        df.iloc[idxs[f"index{suffixes[1]}"].values]
        .rename(columns=lambda x: x + suffixes[1])
        .reset_index(drop=True)
    )

    return pd.concat([left_ivals, right_ivals], axis=1)
