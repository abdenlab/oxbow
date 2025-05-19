from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import Any, Callable, Generator, Iterable

import pyarrow as pa

from oxbow._filetypes import FILETYPE_BY_NAME, FileType
from oxbow._pyarrow import (
    DEFAULT_BATCH_SIZE,
    BatchReaderDataset,
    BatchReaderFragment,
)


class DataSourceMeta(ABCMeta):
    _registry: dict[FileType, type] = {}

    # Invoked when an attribute is not found on the class object
    def __getattr__(cls, key):
        if key.startswith("from_"):
            try:
                file_type_name = key.rsplit("_", 1)[-1].lower()
                file_type = FILETYPE_BY_NAME[file_type_name]
                return {k: v for k, v in cls._registry.items() if issubclass(v, cls)}[
                    file_type
                ]
            except Exception as e:
                raise KeyError(key) from e
        raise AttributeError(f"'{cls.__name__}' object has no attribute '{key}'")


class DataSource(metaclass=DataSourceMeta):
    _schema = None

    @property
    @abstractmethod
    def _scan_kwargs(self) -> dict[str, Any]: ...

    @property
    @abstractmethod
    def _scanner_kwargs(self) -> dict[str, Any]: ...

    @property
    @abstractmethod
    def _schema_kwargs(self) -> dict[str, Any]: ...

    @property
    @abstractmethod
    def _batch_readers(
        self,
    ) -> Iterable[Callable[[list[str] | None, int], pa.RecordBatchReader]]: ...

    def __init_subclass__(cls, file_type: FileType | None = None) -> None:
        if file_type:
            assert file_type not in cls._registry
            cls._scanner_type = file_type.value
            cls._registry[file_type] = cls
        return super().__init_subclass__()

    def __init__(self, uri, opener=None, fields=None):
        self._uri = uri
        self._opener = opener
        self._fields = fields
        super().__init__()

    @property
    def _source(self):
        if self._opener is None:
            return self._uri
        else:
            return self._opener(self._uri)

    @property
    def _scanner(self):
        return self._scanner_type(self._source, **self._scanner_kwargs)

    @property
    def schema(self) -> pa.Schema:
        if not self._schema:
            self._schema = pa.schema(self._scanner.schema(**self._schema_kwargs))
        return self._schema

    @property
    def fields(self) -> dict[str, pa.DataType]:
        return dict(zip(self.schema.names, self.schema.types))

    def select(
        self,
        *regions: str | list[str],
        batch_size: int = DEFAULT_BATCH_SIZE,
        **kwargs,
    ) -> BatchReaderDataset:
        """
        Select a subset of the data source.

        This method creates a new instance of the data source with the same
        parameters, applies any overrides specified by keyword arguments,
        and returns it as a dataset.

        Parameters
        ----------
        *regions
            The regions to select from the data source. This can be a single
            region or a list of regions.
        batch_size
            The size of each batch to be read. Default is DEFAULT_BATCH_SIZE.
        **kwargs
            Additional keyword arguments to override the current instance's
            parameters.

        Returns
        -------
        BatchReaderDataset
            A dataset containing the selected subset of the data source.
        """
        return type(self)(
            self._uri,
            self._opener,
            regions=regions,
            **{
                **self._scan_kwargs,
                **self._scanner_kwargs,
                **self._schema_kwargs,
                **kwargs,
            },
        ).dataset(batch_size=batch_size)

    def batches(
        self, *, batch_size: int = DEFAULT_BATCH_SIZE
    ) -> Generator[pa.RecordBatch]:
        """
        Generate record batches from the data source.

        This method iterates over the data source and yields record batches
        of the specified size.

        Parameters
        ----------
        batch_size
            The size of each batch to be read. Default is DEFAULT_BATCH_SIZE.

        Yields
        ------
        pa.RecordBatch
            A record batch from the data source.
        """
        for reader_factory in self._batch_readers:
            reader = reader_factory(self._fields, batch_size)
            while True:
                try:
                    yield reader.read_next_batch()
                except StopIteration:
                    break

    def fragments(
        self, *, batch_size: int = DEFAULT_BATCH_SIZE
    ) -> list[BatchReaderFragment]:
        """
        Get fragments of the data source.

        Fragments represent parts of the data source that can be processed
        independently. Each fragment is associated with a schema and batch size.

        Parameters
        ----------
        batch_size
            The size of each batch to be read. Default is DEFAULT_BATCH_SIZE.

        Returns
        -------
        list of BatchReaderFragment
            A list of fragments representing parts of the data source.
        """
        return [
            BatchReaderFragment(r, self.schema, batch_size=batch_size)
            for r in self._batch_readers
        ]

    def dataset(self, *, batch_size=DEFAULT_BATCH_SIZE) -> BatchReaderDataset:
        """
        Convert the data source into a dataset.

        A dataset is a collection of fragments that can be processed
        as a single logical entity.

        Parameters
        ----------
        batch_size
            The size of each batch to be read. Default is DEFAULT_BATCH_SIZE.

        Returns
        -------
        BatchReaderDataset
            A dataset representation of the data source.
        """
        return BatchReaderDataset(self.fragments(batch_size=batch_size))

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
            return pl.from_arrow(self.dataset().iter_batches())

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

        # TODO: A less hacky way to generate a "meta" pandas DataFrame for dask
        # without having to materialize a full fragment.
        meta = next(self._fragments[0]._make_batchreader(batch_size=1)).to_pandas()

        # This does a full pass scan over all record batches to find the row
        # offsets of each fragment. While a costly first step, it will endow
        # the Dask DataFrame with "known divisions" which can be exploited for
        # more efficient computations.
        if find_divisions:
            fragment_lengths = [frag.count_rows() for frag in self._fragments]
            row_offsets = [0, *fragment_lengths[:-1]]
            return dd.from_map(
                create_partition,
                self._fragments,
                row_offsets,
                divisions=fragment_lengths,
                meta=meta,
            )

        return dd.from_map(create_partition, self._fragments, meta=meta)

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
        Serialize the dataset as Arrow IPC.

        Returns
        -------
        bytes
            The serialized dataset in Arrow IPC format.
        """
        s = pa.BufferOutputStream()
        with pa.ipc.new_stream(s, self.schema) as writer:
            writer.write_table(self.dataset().to_table())
        buffer = s.getvalue()
        return buffer.to_pybytes()
