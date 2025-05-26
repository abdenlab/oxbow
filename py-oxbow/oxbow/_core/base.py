from __future__ import annotations

from abc import abstractmethod
from typing import Any, Callable, Generator, Iterable, IO, Self
import pathlib

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
        Additional keyword arguments for building the scanner.
    _schema_kwargs : dict
        Additional keyword arguments for assembling the schema.
    _src : Callable[[], IO[Any]]
        A callable that returns the data source.
    _index_src : Callable[[], IO[Any]]
        A callable that returns the index source.
    _batch_size : int
        The size of the batches to be read from the data source.
    _batchreader_builders : Iterable[Callable]
        Callables that create a RecordBatch iterator for each fragment of the
        data source.
    """

    def __init__(
        self,
        source: str | pathlib.Path | Callable[[], IO[Any]],
        index: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        if isinstance(source, (str, pathlib.Path)):
            source = str(source)
            self._src = lambda: source
        elif callable(source):
            self._src = source
        else:
            raise TypeError(
                "`source` must be a str, pathlib.Path, or a callable returning "
                "an IO stream"
            )
        if isinstance(index, (str, pathlib.Path)):
            index = str(index)
            self._index_src = lambda: index
        elif callable(index) or index is None:
            self._index_src = index
        else:
            raise TypeError(
                "`index` must be a str, pathlib.Path, or a callable returning "
                "an IO stream"
            )
        self._batch_size = batch_size

    @property
    def _source(self) -> str | IO[Any]:
        return self._src()

    @property
    def _index(self) -> str | IO[Any]:
        return self._index_src() if self._index_src else None

    @property
    def schema(self) -> pa.Schema:
        return pa.schema(self.scanner().schema(**self._schema_kwargs))

    @property
    def columns(self) -> list[str]:
        return self.schema.names

    @property
    @abstractmethod
    def _batchreader_builders(
        self,
    ) -> Iterable[Callable[[list[str] | None, int], pa.RecordBatchReader]]:
        """
        Callables that generate RecordBatch iterators from the data source.

        Each callable corresponds to a specific fragment of the data source,
        takes column projection and batch size arguments, and returns a
        stream of record batches that cover the fragment.
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
        """
        ...

    def scanner(self) -> Any:
        """
        Create a low-level scanner for the data source.
        """
        return self._scanner_type(self._source, **self._scanner_kwargs)

    def batches(self) -> Generator[pa.RecordBatch]:
        """
        Generate record batches from the data source.

        Yields
        ------
        pa.RecordBatch
            A record batch from the data source.
        """
        for builder in self._batchreader_builders:
            reader = builder(self.columns, self._batch_size)
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
        return [
            BatchReaderFragment(builder, self.schema, batch_size=self._batch_size)
            for builder in self._batchreader_builders
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
            return pl.scan_pyarrow_dataset(self.dataset())
        else:
            batches = list(self.batches())
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
