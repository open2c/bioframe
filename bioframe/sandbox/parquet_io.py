import pandas as pd


def to_parquet(
    pieces,
    outpath,
    row_group_size=None,
    compression="snappy",
    use_dictionary=True,
    version=2.0,
    **kwargs,
):
    """
    Save an iterable of dataframe chunks to a single Apache Parquet file. For
    more info about Parquet, see https://arrow.apache.org/docs/python/parquet.html.

    Parameters
    ----------
    pieces : DataFrame or iterable of DataFrame
        Chunks to write
    outpath : str
        Path to output file
    row_group_size : int
        Number of rows per row group
    compression : {'snappy', 'gzip', 'brotli', 'none'}, optional
        Compression algorithm. Can be set on a per-column basis with a
        dictionary of column names to compression lib.
    use_dictionary : bool, optional
        Use dictionary encoding. Can be set on a per-column basis with a list
        of column names.

    See also
    --------
    pyarrow.parquet.write_table
    pyarrow.parquet.ParquetFile
    fastparquet

    """
    try:
        import pyarrow as pa
        import pyarrow.parquet
    except ImportError:
        raise ImportError("Saving to parquet requires the `pyarrow` package")

    if isinstance(pieces, pd.DataFrame):
        pieces = (pieces,)

    try:
        for i, piece in enumerate(pieces):
            table = pa.Table.from_pandas(piece, preserve_index=False)
            if i == 0:
                writer = pa.parquet.ParquetWriter(
                    outpath,
                    table.schema,
                    compression=compression,
                    use_dictionary=use_dictionary,
                    version=version,
                    **kwargs,
                )
            writer.write_table(table, row_group_size=row_group_size)
    finally:
        writer.close()


def read_parquet(filepath, columns=None, iterator=False, **kwargs):
    """
    Load DataFrames from Parquet files, optionally in pieces.

    Parameters
    ----------
    filepath : str, pathlib.Path, pyarrow.NativeFile, or file-like object
        Readable source. For passing bytes or buffer-like file containing a
        Parquet file, use pyarorw.BufferReader
    columns: list
        If not None, only these columns will be read from the row groups. A
        column name may be a prefix of a nested field, e.g. 'a' will select
        'a.b', 'a.c', and 'a.d.e'
    iterator : boolean, default False
        Return an iterator object that yields row group DataFrames and
        provides the ParquetFile interface.
    use_threads : boolean, default True
        Perform multi-threaded column reads
    memory_map : boolean, default True
        If the source is a file path, use a memory map to read file, which can
        improve performance in some environments

    Returns
    -------
    DataFrame or ParquetFileIterator

    """
    use_threads = kwargs.pop("use_threads", True)

    if not iterator:
        return pd.read_parquet(
            filepath, columns=columns, use_threads=use_threads, **kwargs
        )
    else:
        try:
            from pyarrow.parquet import ParquetFile
        except ImportError:
            raise ImportError(
                "Iterating over Parquet data requires the `pyarrow` package."
            )

        class ParquetFileIterator(ParquetFile):
            def __iter__(self):
                return self

            def __next__(self):
                if not hasattr(self, "_rgid"):
                    self._rgid = 0
                if self._rgid < self.num_row_groups:
                    rg = self.read_row_group(
                        self._rgid,
                        columns=columns,
                        use_threads=use_threads,
                        use_pandas_metadata=True,
                    )
                    self._rgid += 1
                else:
                    raise StopIteration
                return rg.to_pandas()

        return ParquetFileIterator(filepath, **kwargs)
