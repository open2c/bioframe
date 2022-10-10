from collections import OrderedDict
from contextlib import closing
import subprocess
import tempfile
import json
import io

import os
import shutil

import numpy as np
import pandas as pd

try:
    import bbi
except ImportError:
    bbi = None

try:
    import pyBigWig
except ImportError:
    pyBigWig = None

from ..core.stringops import parse_region
from ..core.arrops import argnatsort
from .schemas import SCHEMAS, BAM_FIELDS


__all__ = [
    "read_table",
    "read_chromsizes",
    "read_tabix",
    "read_pairix",
    "read_bam",
    "load_fasta",
    "read_bigwig",
    "to_bigwig",
    "read_bigbed",
    "to_bigbed",
]


def read_table(filepath_or, schema=None, schema_is_strict=False, **kwargs):
    """
    Read a tab-delimited file into a data frame.

    Equivalent to :func:`pandas.read_table` but supports an additional
    `schema` argument to populate column names for common genomic formats.

    Parameters
    ----------
    filepath_or : str, path object or file-like object
        Any valid string path is acceptable. The string could be a URL
    schema : str
        Schema to use for table column names.
    schema_is_strict : bool
        Whether to check if columns are filled with NAs.

    Returns
    -------
    df : pandas.DataFrame of intervals

    """
    kwargs.setdefault("sep", "\t")
    kwargs.setdefault("header", None)
    kwargs.setdefault("index_col", False)
    if isinstance(filepath_or, str) and filepath_or.endswith(".gz"):
        kwargs.setdefault("compression", "gzip")
    if schema is not None:
        try:
            kwargs.setdefault("names", SCHEMAS[schema])
        except (KeyError, TypeError):
            if isinstance(schema, str):
                raise ValueError("TSV schema not found: '{}'".format(schema))
            kwargs.setdefault("names", schema)
    df = pd.read_csv(filepath_or, **kwargs)
    if schema_is_strict:
        if (df.notna().sum(axis=0) == 0).any():
            raise ValueError(
                "one or more columns are all NA,"
                + " check agreement between number of fields in schema"
                + " and number of columns in input file"
            )
    return df


def read_chromsizes(
    filepath_or,
    filter_chroms=True,
    chrom_patterns=(r"^chr[0-9]+$", r"^chr[XY]$", r"^chrM$"),
    natsort=True,
    as_bed=False,
    **kwargs
):
    """
    Read a ``<db>.chrom.sizes`` or ``<db>.chromInfo.txt`` file from the UCSC
    database, where ``db`` is a genome assembly name, as a `pandas.Series`.

    Parameters
    ----------
    filepath_or : str or file-like
        Path or url to text file, or buffer.
    filter_chroms : bool, optional
        Filter for chromosome names given in ``chrom_patterns``.
    chrom_patterns : sequence, optional
        Sequence of regular expressions to capture desired sequence names.
    natsort : bool, optional
        Sort each captured group of names in natural order. Default is True.
    as_bed : bool, optional
        If True, return chromsizes as an interval dataframe (chrom, start, end).
    **kwargs :
        Passed to :func:`pandas.read_csv`

    Returns
    -------
    Series of integer bp lengths indexed by sequence name or an interval dataframe.

    Notes
    -----
    Mention name patterns

    See also
    --------
    * UCSC assembly terminology: <http://genome.ucsc.edu/FAQ/FAQdownloads.html#download9>
    * NCBI assembly terminology: <https://www.ncbi.nlm.nih.gov/grc/help/definitions>

    """
    if isinstance(filepath_or, str) and filepath_or.endswith(".gz"):
        kwargs.setdefault("compression", "gzip")

    chromtable = pd.read_csv(
        filepath_or,
        sep="\t",
        usecols=[0, 1],
        names=["name", "length"],
        dtype={"name": str},
        **kwargs
    )

    if filter_chroms:
        parts = []
        for pattern in chrom_patterns:
            if not len(pattern):
                continue
            part = chromtable[chromtable["name"].str.contains(pattern)]
            if natsort:
                part = part.iloc[argnatsort(part["name"])]
            parts.append(part)
        chromtable = pd.concat(parts, axis=0)

    if as_bed:
        chromtable["start"] = 0
        chromtable = (
            chromtable[["name", "start", "length"]]
            .rename({"name": "chrom", "length": "end"}, axis="columns")
            .reset_index(drop=True)
        )
    else:
        chromtable.index = chromtable["name"].values
        chromtable = chromtable["length"]
    return chromtable


def read_tabix(fp, chrom=None, start=None, end=None):
    """
    Read a tabix-indexed file into dataFrame.
    """
    import pysam

    with closing(pysam.TabixFile(fp)) as f:
        names = list(f.header) or None
        df = pd.read_csv(
            io.StringIO("\n".join(f.fetch(chrom, start, end))),
            sep="\t",
            header=None,
            names=names,
        )
    return df


def read_pairix(
    fp,
    region1,
    region2=None,
    chromsizes=None,
    columns=None,
    usecols=None,
    dtypes=None,
    **kwargs
):
    """
    Read a pairix-indexed file into DataFrame.
    """
    import pypairix
    import cytoolz as toolz

    if dtypes is None:
        dtypes = {}
    f = pypairix.open(fp, "r")

    header = f.get_header()
    if len(header):
        header_groups = toolz.groupby(lambda x: x.split(":")[0], header)
        if "#chromsize" in header_groups and chromsizes is None:
            items = [line.split()[1:] for line in header_groups["#chromsize"]]
            if len(items) and chromsizes is None:
                names, lengths = zip(*((item[0], int(item[1])) for item in items))
                chromsizes = pd.Series(index=names, data=lengths)
        if "#columns" in header_groups and columns is None:
            columns = header_groups["#columns"][0].split()[1:]

    chrom1, start1, end1 = parse_region(region1, chromsizes)
    if region2 is not None:
        chrom2, start2, end2 = parse_region(region2, chromsizes)
    else:
        chrom2, start2, end2 = chrom1, start1, end1

    it = f.query2D(chrom1, start1, end1, chrom2, start2, end2)
    if usecols is not None:
        argusecols = [columns.index(col) for col in usecols]
        records = [(record[i] for i in argusecols) for record in it]
        columns = usecols
    else:
        records = it

    df = pd.DataFrame.from_records(records, columns=columns)
    if columns is not None:
        for col in columns:
            if col in dtypes:
                df[col] = df[col].astype(dtypes[col])
            else:
                df[col] = pd.to_numeric(df[col], "ignore")
    return df


def read_bam(fp, chrom=None, start=None, end=None):
    """
    Read bam records into a DataFrame.
    """
    import pysam

    with closing(pysam.AlignmentFile(fp, "rb")) as f:
        bam_iter = f.fetch(chrom, start, end)
        records = [
            (
                s.qname,
                s.flag,
                s.rname,
                s.pos,
                s.mapq,
                s.cigarstring if s.mapq != 0 else np.nan,
                s.rnext,
                s.pnext,
                s.tlen,
                s.seq,
                s.qual,
                json.dumps(OrderedDict(s.tags)),
            )
            for s in bam_iter
        ]
        df = pd.DataFrame(records, columns=BAM_FIELDS)
    return df


def extract_centromeres(df, schema=None, merge=True):
    """
    Attempts to extract centromere locations from a variety of file formats,
    returning 'chrom', 'start', 'end', 'mid' in a pandas.DataFrame.
    """

    if schema == "centromeres":
        cens = df
    elif schema == "cytoband":
        cens = df[df["gieStain"] == "acen"]
    elif schema == "gap":
        cens = df[df["type"] == "centromere"]
    else:
        raise ValueError('`schema` must be one of {"centromeres", "cytoband", "gap"}.')

    if merge:
        cens = cens.groupby("chrom").agg({"start": np.min, "end": np.max}).reset_index()

    cens["mid"] = (cens["start"] + cens["end"]) // 2
    cens = (
        cens[["chrom", "start", "end", "mid"]]
        .sort_values("chrom")
        .reset_index(drop=True)
    )

    return cens


class PysamFastaRecord(object):
    def __init__(self, ff, ref):
        self.ff = ff
        if ref not in ff.references:
            raise KeyError("Reference name '{}' not found in '{}'".format(ref, ff))
        self.ref = ref

    def __getitem__(self, key):
        if isinstance(key, slice):
            start, stop = key.start, key.stop
        else:
            start = key
            stop = key + 1
        return self.ff.fetch(self.ref, start, stop)


def load_fasta(filepath_or, engine="pysam", **kwargs):
    """
    Load lazy fasta sequences from an indexed fasta file (optionally compressed)
    or from a collection of uncompressed fasta files.

    Parameters
    ----------
    filepath_or : str or iterable
        If a string, a filepath to a single `.fa` or `.fa.gz` file. Assumed to
        be accompanied by a `.fai` index file. Depending on the engine, the
        index may be created on the fly, and some compression formats may not
        be supported. If not a string, an iterable of fasta file paths each
        assumed to contain a single sequence.
    engine : {'pysam', 'pyfaidx'}, optional
        Module to use for loading sequences.
    kwargs : optional
        Options to pass to ``pysam.FastaFile`` or ``pyfaidx.Fasta``.

    Returns
    -------
    OrderedDict of (lazy) fasta records.

    Notes
    -----
    * pysam/samtools can read .fai and .gzi indexed files, I think.
    * pyfaidx can handle uncompressed and bgzf compressed files.

    """
    is_multifile = not isinstance(filepath_or, str)
    records = OrderedDict()

    engine = engine.lower()

    if engine == "pysam":
        try:
            import pysam
        except ImportError:
            raise ImportError("pysam is required to use engine='pysam'")

        if is_multifile:
            for onefile in filepath_or:
                ff = pysam.FastaFile(onefile, **kwargs)
                name = ff.references[0]
                records[name] = PysamFastaRecord(ff, name)
        else:
            ff = pysam.FastaFile(filepath_or, **kwargs)
            for name in ff.references:
                records[name] = PysamFastaRecord(ff, name)

    elif engine == "pyfaidx":
        try:
            import pyfaidx
        except ImportError:
            raise ImportError("pyfaidx is required to use engine='pyfaidx'")

        if is_multifile:
            for onefile in filepath_or:
                ff = pyfaidx.Fasta(onefile, **kwargs)
                name = next(iter(ff.keys()))
                records[name] = ff[name]
        else:
            ff = pyfaidx.Fasta(filepath_or, **kwargs)
            for name in ff.keys():
                records[name] = ff[name]

    else:
        raise ValueError("engine must be 'pysam' or 'pyfaidx'")

    return records


def read_bigwig(path, chrom, start=None, end=None, engine="auto"):
    """
    Read intervals from a bigWig file.

    Parameters
    ----------
    path : str
        Path or URL to a bigWig file
    chrom : str
    start, end : int, optional
        Start and end coordinates. Defaults to 0 and chromosome length.
    engine : {"auto", "pybbi", "pybigwig"}
        Library to use for querying the bigWig file.

    Returns
    -------
    DataFrame

    """
    engine = engine.lower()

    if engine == "auto":
        if bbi is None and pyBigWig is None:
            raise ImportError(
                "read_bigwig requires either the pybbi or pyBigWig package"
            )
        elif bbi is not None:
            engine = "pybbi"
        else:
            engine = "pybigwig"

    if engine in ("pybbi", "bbi"):
        if start is None:
            start = 0
        if end is None:
            end = -1
        with bbi.open(path) as f:
            df = f.fetch_intervals(chrom, start=start, end=end)

    elif engine == "pybigwig":
        f = pyBigWig.open(path)
        if start is None:
            start = 0
        if end is None:
            end = f.chroms()[chrom]
        ivals = f.intervals(chrom, start, end)
        df = pd.DataFrame(ivals, columns=["start", "end", "value"])
        df.insert(0, "chrom", chrom)

    else:
        raise ValueError(
            "engine must be 'auto', 'pybbi' or 'pybigwig'; got {}".format(engine)
        )

    return df


def read_bigbed(path, chrom, start=None, end=None, engine="auto"):
    """
    Read intervals from a bigBed file.

    Parameters
    ----------
    path : str
        Path or URL to a bigBed file
    chrom : str
    start, end : int, optional
        Start and end coordinates. Defaults to 0 and chromosome length.
    engine : {"auto", "pybbi", "pybigwig"}
        Library to use for querying the bigBed file.

    Returns
    -------
    DataFrame

    """
    engine = engine.lower()

    if engine == "auto":
        if bbi is None and pyBigWig is None:
            raise ImportError(
                "read_bigbed requires either the pybbi or pyBigWig package"
            )
        elif bbi is not None:
            engine = "pybbi"
        else:
            engine = "pybigwig"

    if engine in ("pybbi", "bbi"):
        if start is None:
            start = 0
        if end is None:
            end = -1
        with bbi.open(path) as f:
            df = f.fetch_intervals(chrom, start=start, end=end)

    elif engine == "pybigwig":
        f = pyBigWig.open(path)
        if start is None:
            start = 0
        if end is None:
            end = f.chroms()[chrom]
        ivals = f.entries(chrom, start, end)
        df = pd.DataFrame(ivals, columns=["start", "end", "rest"])
        df.insert(0, "chrom", chrom)

    else:
        raise ValueError(
            "engine must be 'auto', 'pybbi' or 'pybigwig'; got {}".format(engine)
        )

    return df


def to_bigwig(df, chromsizes, outpath, value_field=None, path_to_binary=None):
    """
    Save a bedGraph-like dataframe as a binary BigWig track.

    Parameters
    ----------
    df : pandas.DataFrame
        Data frame with columns 'chrom', 'start', 'end' and one or more value
        columns
    chromsizes : pandas.Series
        Series indexed by chromosome name mapping to their lengths in bp
    outpath : str
        The output BigWig file path
    value_field : str, optional
        Select the column label of the data frame to generate the track. Default
        is to use the fourth column.
    path_to_binary : str, optional
        Provide system path to the bedGraphToBigWig binary.

    """

    if path_to_binary is None:
        cmd = "bedGraphToBigWig"
        try:
            assert shutil.which(cmd) is not None
        except Exception as e:
            raise ValueError(
                "bedGraphToBigWig is not present in the current environment. "
                "Pass it as 'path_to_binary' parameter to bioframe.to_bigwig or "
                "install it with, for example, conda install -y -c bioconda ucsc-bedgraphtobigwig "
            )
    elif path_to_binary.endswith("bedGraphToBigWig"):
        if not os.path.isfile(path_to_binary) and os.access(path_to_binary, os.X_OK):
            raise ValueError(
                f"bedGraphToBigWig is absent in the provided path or cannot be executed: {path_to_binary}. "
            )
        cmd = path_to_binary
    else:
        cmd = os.path.join(path_to_binary, "bedGraphToBigWig")
        if not os.path.isfile(cmd) and os.access(cmd, os.X_OK):
            raise ValueError(
                f"bedGraphToBigWig is absent in the provided path or cannot be executed: {path_to_binary}. "
            )

    is_bedgraph = True
    for col in ["chrom", "start", "end"]:
        if col not in df.columns:
            is_bedgraph = False
    if len(df.columns) < 4:
        is_bedgraph = False

    if not is_bedgraph:
        raise ValueError(
            "A bedGraph-like DataFrame is required, got {}".format(df.columns)
        )

    if value_field is None:
        value_field = df.columns[3]

    columns = ["chrom", "start", "end", value_field]
    bg = df[columns].copy()
    bg["chrom"] = bg["chrom"].astype(str)
    bg = bg.sort_values(["chrom", "start", "end"])

    with tempfile.NamedTemporaryFile(suffix=".bg") as f, tempfile.NamedTemporaryFile(
       "wt", suffix=".chrom.sizes"
    ) as cs:

        chromsizes.to_csv(cs, sep="\t", header=False)
        cs.flush()

        bg.to_csv(
            f.name, sep="\t", columns=columns, index=False, header=False, na_rep="nan"
        )

        p = subprocess.run(
            [cmd, f.name, cs.name, outpath],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    return p


def to_bigbed(df, chromsizes, outpath, schema="bed6", path_to_binary=None):
    """
    Save a bedGraph-like dataframe as a binary BigWig track.

    Parameters
    ----------
    df : pandas.DataFrame
        Data frame with columns 'chrom', 'start', 'end' and one or more value
        columns
    chromsizes : pandas.Series
        Series indexed by chromosome name mapping to their lengths in bp
    outpath : str
        The output BigWig file path
    value_field : str, optional
        Select the column label of the data frame to generate the track. Default
        is to use the fourth column.
    path_to_binary : str, optional
        Provide system path to the bedGraphToBigWig binary.

    """

    if path_to_binary is None:
        cmd = "bedToBigBed"
        try:
            assert shutil.which(cmd) is not None
        except Exception as e:
            raise ValueError(
                "bedToBigBed is not present in the current environment. "
                "Pass it as 'path_to_binary' parameter to bioframe.to_bigbed or "
                "install it with, for example, conda install -y -c bioconda ucsc-bedtobigbed "
            )
    elif path_to_binary.endswith("bedToBigBed"):
        if not os.path.isfile(path_to_binary) and os.access(path_to_binary, os.X_OK):
            raise ValueError(
                f"bedToBigBed is absent in the provided path or cannot be executed: {path_to_binary}. "
            )
        cmd = path_to_binary
    else:
        cmd = os.path.join(path_to_binary, "bedGraphToBigWig")
        if not os.path.isfile(cmd) and os.access(cmd, os.X_OK):
            raise ValueError(
                f"bedToBigBed is absent in the provided path  or cannot be executed: {path_to_binary}. "
            )

    is_bed6 = True
    for col in ["chrom", "start", "end", "name", "score", "strand"]:
        if col not in df.columns:
            is_bed6 = False
    if len(df.columns) < 6:
        is_bed6 = False

    if not is_bed6:
        raise ValueError("A bed6-like DataFrame is required, got {}".format(df.columns))

    columns = ["chrom", "start", "end", "name", "score", "strand"]
    bed = df[columns].copy()
    bed["chrom"] = bed["chrom"].astype(str)
    bed = bed.sort_values(["chrom", "start", "end"])

    with tempfile.NamedTemporaryFile(suffix=".bed") as f, tempfile.NamedTemporaryFile(
        "wt", suffix=".chrom.sizes"
    ) as cs:

        chromsizes.to_csv(cs, sep="\t", header=False)
        cs.flush()

        bed.to_csv(
            f.name, sep="\t", columns=columns, index=False, header=False, na_rep="nan"
        )

        p = subprocess.run(
            [cmd, "-type={}".format(schema), f.name, cs.name, outpath],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    return p
