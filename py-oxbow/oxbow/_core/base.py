from __future__ import annotations

import pathlib
import warnings
from abc import abstractmethod
from typing import IO, Any, Callable, Generator, Literal

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

from urllib.parse import urlparse

import fsspec
import pyarrow as pa

from oxbow._pyarrow import (
    DEFAULT_BATCH_SIZE,
    BatchReaderDataset,
    BatchReaderFragment,
)


class DataSource:
    """
    Base class for data sources.

    Attributes
    ----------
    _scanner_type : type
        The scanner type used for reading the data source.
    _scanner_kwargs : dict
        Keyword arguments for building the scanner (includes schema params).
    _src : Callable[[], IO[bytes]]
        A callable that returns the data source.
    _index_src : Callable[[], IO[bytes]]
        A callable that returns the index source.
    _batch_size : int
        The size of the batches to be read from the data source.
    """

    _scanner_type: type
    _scanner_kwargs: dict[str, Any] = {}

    def __init__(
        self,
        source: str | Callable[[], IO[bytes] | str],
        index: str | Callable[[], IO[bytes] | str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        self._src = source
        self._index_src = index
        self._batch_size = batch_size

    @property
    def _source(self) -> IO[bytes] | str:
        return self._src() if callable(self._src) else self._src

    @property
    def _index(self) -> IO[bytes] | str | None:
        return (
            self._index_src()
            if (self._index_src and callable(self._index_src))
            else self._index_src
        )

    @property
    def schema(self) -> pa.Schema:
        """The arrow schema of the projection."""
        return pa.schema(self.scanner().schema())

    @property
    def columns(self) -> list[str]:
        """The top-level column names of the projection."""
        return self.schema.names

    def _make_reader(self, columns, batch_size, region=None):
        """Return a RecordBatchReader for a full scan or regional query.

        The default implementation delegates to ``scanner.scan()`` for full
        scans and to :meth:`_scan_query` for regional queries.  Subclasses
        may override this entirely when the scanner API differs from the
        common ``columns`` / ``batch_size`` convention.

        Parameters
        ----------
        columns : list[str] | None
            Column projection, or None for all columns.
        batch_size : int
            Number of records per batch.
        region : str, optional [default: None]
            A genomic region string, e.g. ``"chr1:1000-2000"``.

        Returns
        -------
        RecordBatchReader
            A RecordBatchReader for the specified scan.
        """
        scanner = self.scanner()
        if region is not None:
            return self._scan_query(scanner, region, columns, batch_size)
        return scanner.scan(columns=columns, batch_size=batch_size)

    @abstractmethod
    def _scan_query(self, scanner, region, columns, batch_size):
        """Return a RecordBatchReader for a range query.

        Subclasses implement this to handle format-specific index and query
        semantics (e.g. passing an index object, handling unmapped reads).

        Parameters
        ----------
        scanner
            The low-level scanner instance (already constructed).
        region : str
            A genomic region string, e.g. ``"chr1:1000-2000"``.
        columns : list[str] | None
            Column projection, or None for all columns.
        batch_size : int
            Number of records per batch.

        Returns
        -------
        RecordBatchReader
            A RecordBatchReader for the specified range query.
        """
        ...

    @abstractmethod
    def regions(self, regions: str | list[str]) -> Self:
        """
        Query one or more genomic ranges within the data source.

        This method creates a new instance of the data source with the same
        parameters, overriding the regions to select from the data source.

        Parameters
        ----------
        regions: str | list[str]
            The regions to select from the data source. This can be a single
            region or a list of regions.

        Returns
        -------
        DataSource

        Notes
        -----
        Genomic range strings can be in the following formats:

        - UCSC-style ``"chr:start-end"``: intepreted using the coordinate
          system of the data source.
        - Bracket-style ``"chr:[start,end]"``: explicitly 1-based, end-inclusive.
        - Bracket-style ``"chr:[start,end)"``: explicitly 0-based, end-exclusive.
        """
        ...

    def scanner(self) -> Any:
        """
        Create a low-level scanner for the data source.
        """
        return self._scanner_type(self._source, **self._scanner_kwargs)

    def batches(self) -> Generator:
        """
        Generate record batches from the data source.

        Yields
        ------
        RecordBatch
            A record batch from the data source.
        """
        for region in self._regions or [None]:
            reader = self._make_reader(self.columns, self._batch_size, region)
            while True:
                try:
                    yield reader.read_next_batch()
                except StopIteration:
                    break

    def fragments(self) -> list[BatchReaderFragment]:
        """
        Get fragments of the data source.

        Fragments represent parts of the data source that can be processed
        independently.

        Returns
        -------
        list of BatchReaderFragment
            A list of fragments representing parts of the data source.
        """
        schema = self.schema
        regions = self._regions or [None]
        return [
            BatchReaderFragment(
                lambda columns, batch_size, r=region: self._make_reader(
                    columns, batch_size, r
                ),
                schema,
                batch_size=self._batch_size,
                tokenize=(
                    # deterministic only if source and index are str
                    self._source,
                    self._index,
                    region,
                    self._scanner_kwargs,
                    self._batch_size,
                ),
            )
            for region in regions
        ]

    def dataset(self) -> BatchReaderDataset:
        """
        Convert the data source into a dataset.

        A dataset is a collection of fragments that can be processed
        as a single logical entity.

        Returns
        -------
        BatchReaderDataset
            A dataset representation of the data source.
        """
        return BatchReaderDataset(self.fragments())

    def to_pandas(self):
        """
        Convert the dataset to a Pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame representation of the dataset.
        """
        return self.dataset().to_table().to_pandas()

    pd = to_pandas

    def to_polars(self, lazy=False):
        """
        Convert the data source to a Polars DataFrame or LazyFrame.

        Parameters
        ----------
        lazy : bool, optional [default: False]
            If True, returns a LazyFrame.

        Returns
        -------
        polars.DataFrame | polars.LazyFrame
            A polars representation of the data source.
        """
        import polars as pl

        if lazy:
            from polars.io.plugins import register_io_source

            regions = self._regions or [None]
            arrow_schema = self.schema
            default_batch_size = self._batch_size
            polars_schema = pl.from_arrow(arrow_schema.empty_table()).schema

            def io_source(with_columns, predicate, n_rows, batch_size):
                rows_read = 0
                bs = batch_size or default_batch_size
                for region in regions:
                    reader = self._make_reader(with_columns, bs, region)
                    while True:
                        try:
                            batch = reader.read_next_batch()
                        except StopIteration:
                            break
                        df = pl.from_arrow(batch)
                        if predicate is not None:
                            df = df.filter(predicate)
                        if n_rows is not None:
                            remaining = n_rows - rows_read
                            if remaining <= 0:
                                return
                            if len(df) > remaining:
                                df = df.head(remaining)
                        rows_read += len(df)
                        yield df
                        if n_rows is not None and rows_read >= n_rows:
                            return

            return register_io_source(io_source=io_source, schema=polars_schema)
        else:
            batches = [pa.record_batch(b) for b in self.batches()]
            if not batches:
                return pl.from_arrow(self.schema.empty_table())
            return pl.from_arrow(batches)

    pl = to_polars

    def to_dask(self, find_divisions=False):
        """
        Convert the data source to a Dask DataFrame.

        Parameters
        ----------
        find_divisions : bool, optional
            If True, find divisions for the Dask DataFrame, by default False.

        Returns
        -------
        dask.dataframe.DataFrame
            A Dask DataFrame representation of the data source.
        """
        import dask.dataframe as dd
        import pandas as pd

        def create_partition(fragment, row_offset=None, columns=None):
            df = fragment.to_table(columns=columns).to_pandas()
            if row_offset is not None:
                df = df.set_index(pd.RangeIndex(row_offset, row_offset + len(df)))
            return df

        # Generate a "meta" pandas DataFrame which serves as a schema for dask
        meta = self.schema.empty_table().to_pandas()

        # This does a full-pass scan over all record batches to find the row
        # offsets of each fragment. While a costly first step, it will endow
        # the Dask DataFrame with "known divisions" which can be exploited for
        # more efficient computations.
        fragments = self.fragments()
        if find_divisions:
            fragment_lengths = [frag.count_rows() for frag in fragments]
            row_offsets = [0, *fragment_lengths[:-1]]
            return dd.from_map(
                create_partition,
                fragments,
                row_offsets,
                divisions=fragment_lengths,
                meta=meta,
            )
        else:
            return dd.from_map(create_partition, fragments, meta=meta)

    dd = to_dask

    def to_duckdb(self, conn):
        """
        Convert the data source into a DuckDB Relation.

        Parameters
        ----------
        conn : duckdb.DuckDBPyConnection
            The DuckDB connection.

        Returns
        -------
        duckdb.DuckDBPyRelation
            A DuckDB Relation representation of the data source.
        """
        return conn.from_arrow(self.dataset())

    def to_ipc(self) -> bytes:
        """
        Serialize the data source as an Arrow IPC stream.

        Returns
        -------
        bytes
            The serialized data source in Arrow IPC format.
        """
        s = pa.BufferOutputStream()
        with pa.ipc.new_stream(s, self.schema) as writer:
            writer.write_table(self.dataset().to_table())
        buffer = s.getvalue()
        return buffer.to_pybytes()


def prepare_source_and_index(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    compression: Literal["infer", "bgzf", "gzip", None] = "infer",
) -> tuple[str | Callable[[], IO[bytes]], str | Callable[[], IO[bytes]] | None, bool]:
    if isinstance(source, (str, pathlib.Path)):
        source = str(source)
        use_fsspec = (scheme := urlparse(source).scheme) and scheme in (
            "http",
            "https",
            "ftp",
            "s3",
            "file",
        )
        if compression == "infer":
            bgzf_compressed = str(source).endswith(".gz") or str(source).endswith(
                ".bgz"
            )
        elif compression == "bgzf":
            bgzf_compressed = True
        elif compression == "gzip":
            bgzf_compressed = False
            use_fsspec = True
        else:
            bgzf_compressed = False

        if use_fsspec:
            src = lambda: fsspec.open(  # noqa: E731
                source,
                mode="rb",
                compression=compression if compression == "gzip" else None,
            ).open()
        else:
            src = source
    elif callable(source):
        src = source
        if compression == "infer":
            warnings.warn(
                "Compression inference is not supported for callable sources. "
                "Assuming bytestream returned by source is uncompressed."
            )
            bgzf_compressed = False
        elif compression == "bgzf":
            bgzf_compressed = True
        elif compression == "gzip":
            raise ValueError(
                "'gzip' compression is not supported for callable sources. "
                "The callable should handle decompression in this case."
            )
        else:
            bgzf_compressed = False
    else:
        raise TypeError(
            "`source` must be a str, pathlib.Path, or a callable returning "
            "an IO byte stream"
        )

    if isinstance(index, (str, pathlib.Path)):
        index = str(index)
        if (scheme := urlparse(index).scheme) and scheme in (
            "http",
            "https",
            "ftp",
            "s3",
            "file",
        ):
            index_src = lambda: fsspec.open(index, mode="rb").open()  # noqa: E731
        else:
            index_src = index
    elif callable(index) or index is None:
        index_src = index
    else:
        raise TypeError(
            "`index` must be a str, pathlib.Path, or a callable returning "
            "an IO byte stream"
        )

    return src, index_src, bgzf_compressed
