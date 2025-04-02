from __future__ import annotations

from itertools import chain
from functools import partial
from typing import Iterator, Callable, Final

import pyarrow as pa
import pyarrow.dataset as ds
from pyarrow.dataset import Fragment, Dataset, Scanner

DEFAULT_BATCH_SIZE: Final = 2**17
DEFAULT_BATCH_READAHEAD: Final = 16
DEFAULT_FRAGMENT_READAHEAD: Final = 4

RecordBatchIter = pa.RecordBatchReader | Iterator[pa.RecordBatch]


class BatchReaderFragment(Fragment):
    """
    A Fragment that emits RecordBatches from a reproducible source.

    To provide stateless replay, a new record batch iterator over the same
    records is constructed whenever a scanner is requested.
    """

    def __init__(
        self,
        make_batchreader: Callable[[list[str] | None, int], RecordBatchIter],
        schema: pa.Schema,
        batch_size: int = DEFAULT_BATCH_SIZE,
        partition_expression: ds.Expression = None,
    ):
        """
        Create a BatchReaderFragment from a BatchReader factory function.

        Parameters
        ----------
        make_batchreader : Callable
            A function that recreates a specific stream of record batches.
        schema : pa.Schema
            The schema of the RecordBatches.
        batch_size : int, optional
            The maximum row count for scanned record batches.
        partition_expression : ds.Expression, optional
            A partition expression for the fragment.

        Notes
        -----
        The `make_batchreader` function should accept the following arguments
        and return a pyarrow.RecordBatchReader:

        * columns: list[str] | None
            The columns to project.
        * batch_size: int
            The maximum number of rows per batch.
        """
        self._make_batchreader = make_batchreader
        self._schema = schema
        self._batch_size = batch_size
        if partition_expression is not None:
            self._partition_expression = partition_expression
        else:
            self._partition_expression = ds.scalar(True)

    @property
    def schema(self) -> pa.Schema:
        """
        Returns
        -------
        schema : pyarrow.Schema
            The schema of the RecordBatches in the fragment.
        """
        return self._schema

    @property
    def partition_expression(self) -> ds.Expression:
        """
        Returns
        -------
        partition_expression : pyarrow.dataset.Expression
            An expression that evaluates to true for all data viewed by this fragment.
        """
        return self._partition_expression

    def scanner(
        self,
        schema: pa.Schema | None = None,
        columns: list[str] | None = None,
        filter: ds.Expression | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
        batch_readahead: int = DEFAULT_BATCH_READAHEAD,
        fragment_readahead: int = DEFAULT_FRAGMENT_READAHEAD,
        fragment_scan_options: ds.FragmentScanOptions | None = None,
        use_threads: bool = True,
        memory_pool: ds.MemoryPool = None,
    ) -> ds.Scanner:
        """
        Build a scan operation against the fragment.

        Parameters
        ----------
        schema : pa.Schema, optional
            The schema to use for scanning. This is used to unify a Fragment to
            its Dataset's schema. If not specified this will use the Fragment's
            physical schema.
        columns : list[str], optional
            Names of columns to project. By default all of the available
            columns are projected.
        filter : pa.Expression, optional
            A filter expression. Scan will return only the rows matching the
            filter.
        batch_size : int, optional
            The maximum row count for scanned record batches.
        batch_readahead : int, default 16
            The number of batches to read ahead in a file.
        fragment_readahead : int, default 4
            The number of fragments/files to read ahead.
        fragment_scan_options : FragmentScanOptions, default None
            Options specific to a particular scan and fragment type, which
            can change between different scans of the same dataset.
        use_threads : bool, default True
            If enabled, then maximum parallelism will be used determined by
            the number of available CPU cores.
        memory_pool : MemoryPool, default None
            For memory allocations, if required. If not specified, uses the
            default pool.

        Returns
        -------
        Scanner

        Notes
        -----
        The list of columns may include special fields:

        * ``__batch_index``: Index of the batch within the fragment.
        * ``__fragment_index``: Index of the fragment within the dataset.
        * ``__last_in_fragment``: Whether the batch is last in fragment.
        * ``__filename``: Name of the source file or a description of the
            source fragment.
        """
        # Apply column projection if specified
        schema = schema or self._schema
        batches = self._make_batchreader(fields=columns, batch_size=batch_size)

        # The scanner constructor treats a RecordBatchReader differently than
        # an opaque iterator, so we wrap it in an iterator for consistency
        def _iterator(batches):
            for batch in batches:
                yield batch

        if isinstance(batches, pa.RecordBatchReader):
            batches = _iterator(batches)

        # Make a Scanner from the batches, applying filter if specified
        return Scanner.from_batches(
            source=batches,
            schema=schema or self.schema,
            columns=columns,
            filter=filter,
            batch_size=batch_size,
            batch_readahead=batch_readahead,
            fragment_readahead=fragment_readahead,
            fragment_scan_options=fragment_scan_options,
            use_threads=use_threads,
            memory_pool=memory_pool,
        )

    def to_batches(
        self,
        schema: pa.Schema | None = None,
        columns: list[str] | None = None,
        filter: ds.Expression | None = None,
        batch_size: int | None = DEFAULT_BATCH_SIZE,
        batch_readahead: int = DEFAULT_BATCH_READAHEAD,
        fragment_readahead: int = DEFAULT_FRAGMENT_READAHEAD,
        fragment_scan_options: ds.FragmentScanOptions | None = None,
        use_threads: bool = True,
        memory_pool: ds.MemoryPool = None,
    ) -> Iterator[pa.RecordBatch]:
        """
        Scan and read the fragment as materialized record batches.

        Projections and filters are applied if specified.

        Parameters
        ----------
        schema : pa.Schema, optional
            The schema to use for scanning. This is used to unify a Fragment to
            its Dataset's schema. If not specified this will use the Fragment's
            physical schema.
        columns : list[str], optional
            Names of columns to project. By default all of the available
            columns are projected.
        filter : pa.Expression, optional
            A filter expression. Scan will return only the rows matching the
            filter.
        batch_size : int, optional
            The maximum row count for scanned record batches.
        batch_readahead : int, default 16
            The number of batches to read ahead in a file.
        fragment_readahead : int, default 4
            The number of fragments/files to read ahead.
        fragment_scan_options : FragmentScanOptions, default None
            Options specific to a particular scan and fragment type, which
            can change between different scans of the same dataset.
        use_threads : bool, default True
            If enabled, then maximum parallelism will be used determined by
            the number of available CPU cores.
        memory_pool : MemoryPool, default None
            For memory allocations, if required. If not specified, uses the
            default pool.

        Returns
        -------
        Iterator[RecordBatch]

        Notes
        -----
        The list of columns may include special fields:

        * ``__batch_index``: Index of the batch within the fragment.
        * ``__fragment_index``: Index of the fragment within the dataset.
        * ``__last_in_fragment``: Whether the batch is last in fragment.
        * ``__filename``: Name of the source file or a description of the
            source fragment.
        """
        return self.scanner(
            schema=schema or self.schema,
            columns=columns,
            filter=filter,
            batch_size=batch_size,
            batch_readahead=batch_readahead,
            fragment_readahead=fragment_readahead,
            fragment_scan_options=fragment_scan_options,
            use_threads=use_threads,
            memory_pool=memory_pool,
        ).to_batches()

    def iter_batches(self, columns=None, batch_size=DEFAULT_BATCH_SIZE):
        return self._make_batchreader(fields=columns, batch_size=batch_size)

    def __dask_tokenize__(self):
        from dask.base import normalize_token

        return (
            normalize_token(self.__class__),
            self._schema,
            self._partition_expression,
            self._batch_size,
        )


class BatchReaderDataset(Dataset):
    """A PyArrow Dataset composed of one or more BatchReaderFragments."""

    def __init__(
        self,
        fragments: list[BatchReaderFragment],
        partition_expression: ds.Expression = None,
    ):
        """
        Create a BatchReaderDataset from a list of BatchReaderFragments.

        Parameters
        ----------
        fragments : list[BatchReaderFragment]
            The list of fragments.

        partition_expression : ds.Expression, optional
            A partition expression for the dataset.
        """
        self._fragments = fragments
        self._schema = self._fragments[0].schema
        self._scan_options = {}
        if partition_expression is not None:
            self._partition_expression = partition_expression
        else:
            self._partition_expression = ds.scalar(True)

    @property
    def partition_expression(self) -> ds.Expression:
        """
        Returns
        -------
        partition_expression : pyarrow.dataset.Expression
            An expression that evaluates to true for all data viewed by this dataset.
        """
        return self._partition_expression

    @property
    def schema(self) -> pa.Schema:
        """
        Returns
        -------
        schema : pyarrow.Schema
            The schema of the RecordBatches in the dataset.
        """
        return self._schema

    def get_fragments(
        self, filter: ds.Expression | None = None
    ) -> Iterator[BatchReaderFragment]:
        """
        Return an iterator over fragments.

        Parameters
        ----------
        filter : pyarrow.dataset.Expression, optional
            An expression to filter fragments. Not yet implemented.

        Returns
        -------
        Iterator[BatchReaderFragment]

        Notes
        -----
        ``filter`` here is meant to be applied at the fragment level via
        comparision with the parition_expression, not at the row level.
        """
        if filter is None:
            for fragment in self._fragments:
                yield fragment
        else:
            # https://github.com/apache/arrow/blob/76fa19e61af25d124ec0af5e543110a4672088db/cpp/src/arrow/dataset/dataset.cc#L204C1-L210C1
            raise NotImplementedError("Fragment-level filtering is not yet implemented")

    def scanner(
        self,
        columns: list[str] | None = None,
        filter: ds.Expression | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
        batch_readahead: int = DEFAULT_BATCH_READAHEAD,
        fragment_readahead: int = DEFAULT_FRAGMENT_READAHEAD,
        fragment_scan_options: ds.FragmentScanOptions | None = None,
        use_threads: bool | None = True,
        memory_pool: ds.MemoryPool | None = None,
    ) -> ds.Scanner:
        """
        Build a scan operation against the dataset.

        This scanner chains the record batches from all fragments together
        and applies column projection and row filtering.

        Parameters
        ----------
        columns : list[str], optional
            Names of columns to project. By default all of the available
            columns are projected.
        filter : pa.Expression, optional
            A filter expression. Scan will return only the rows matching the
            filter.
        batch_size : int, optional
            The maximum row count for scanned record batches.
        batch_readahead : int, default 16
            The number of batches to read ahead in a file.
        fragment_readahead : int, default 4
            The number of fragments/files to read ahead.
        fragment_scan_options : FragmentScanOptions, default None
            Options specific to a particular scan and fragment type, which
            can change between different scans of the same dataset.
        use_threads : bool, default True
            If enabled, then maximum parallelism will be used determined by
            the number of available CPU cores.
        memory_pool : MemoryPool, default None
            For memory allocations, if required. If not specified, uses the
            default pool.

        Returns
        -------
        Scanner

        Notes
        -----
        The list of columns may include special fields:

        * ``__batch_index``: Index of the batch within the fragment.
        * ``__fragment_index``: Index of the fragment within the dataset.
        * ``__last_in_fragment``: Whether the batch is last in fragment.
        * ``__filename``: Name of the source file or a description of the
            source fragment.
        """
        # TODO: Prune fragments using their partition expressions.

        # Chain all the fragments' record batch iterators together; don't
        # apply any filter yet. No batches should get materialized.
        batch_iter = chain.from_iterable(
            fragment.iter_batches(
                columns=columns,
                batch_size=batch_size,
            )
            for fragment in self._fragments
        )

        # Apply the row filter via the scanner.
        return Scanner.from_batches(
            source=batch_iter,
            schema=self.schema,
            columns=columns,
            filter=filter,
            batch_size=batch_size,
            batch_readahead=batch_readahead,
            fragment_readahead=fragment_readahead,
            fragment_scan_options=fragment_scan_options,
            use_threads=use_threads,
            memory_pool=memory_pool,
        )

    def to_batches(
        self,
        columns: list[str] | None = None,
        filter: ds.Expression | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
        batch_readahead: int = DEFAULT_BATCH_READAHEAD,
        fragment_readahead: int = DEFAULT_FRAGMENT_READAHEAD,
        fragment_scan_options: ds.FragmentScanOptions | None = None,
        use_threads: bool | None = True,
        memory_pool: ds.MemoryPool | None = None,
    ) -> Iterator[pa.RecordBatch]:
        """
        Read the dataset as materialized record batches.

        Parameters
        ----------
        columns : list[str], optional
            Names of columns to project. By default all of the available
            columns are projected.
        filter : pa.Expression, optional
            A filter expression. Scan will return only the rows matching the
            filter.
        batch_size : int, optional
            The maximum row count for scanned record batches.
        batch_readahead : int, default 16
            The number of batches to read ahead in a file.
        fragment_readahead : int, default 4
            The number of fragments/files to read ahead.
        fragment_scan_options : FragmentScanOptions, default None
            Options specific to a particular scan and fragment type, which
            can change between different scans of the same dataset.
        use_threads : bool, default True
            If enabled, then maximum parallelism will be used determined by
            the number of available CPU cores.
        memory_pool : MemoryPool, default None
            For memory allocations, if required. If not specified, uses the
            default pool.

        Returns
        -------
        Iterator[RecordBatch]

        Notes
        -----
        The list of columns may include special fields:

        * ``__batch_index``: Index of the batch within the fragment.
        * ``__fragment_index``: Index of the fragment within the dataset.
        * ``__last_in_fragment``: Whether the batch is last in fragment.
        * ``__filename``: Name of the source file or a description of the
            source fragment.
        """
        return self.scanner(
            columns=columns,
            filter=filter,
            batch_size=batch_size,
            batch_readahead=batch_readahead,
            fragment_readahead=fragment_readahead,
            fragment_scan_options=fragment_scan_options,
            use_threads=use_threads,
            memory_pool=memory_pool,
        ).to_batches()

    def iter_batches(self, columns=None, batch_size=DEFAULT_BATCH_SIZE):
        return chain.from_iterable(
            fragment.iter_batches(
                columns=columns,
                batch_size=batch_size,
            )
            for fragment in self._fragments
        )

    def filter(self, expression: ds.Expression) -> "BatchReaderDataset":
        """
        Apply a row-level filter expression and return a filtered dataset.

        Parameters
        ----------
        expression : pa.Expression
            The filter expression.

        Returns
        -------
        filtered_dataset : RecordBatchDataset
            The filtered dataset.

        Notes
        -----
        This performs a row-level filter, not a fragment-level filter.
        The new filter is applied lazily on the underlying record batches by
        using a wrapper factory function to make a new BatchReaderFragment.
        """
        new_filter = expression
        current_filter = self._scan_options.get("filter")
        if current_filter is not None and new_filter is not None:
            new_filter = current_filter & new_filter

        new_fragments = []

        def filter_batches(fragment, columns, batch_size):
            batches = fragment.to_batches(
                schema=self.schema,
                columns=columns,
                filter=new_filter,
                batch_size=batch_size,
            )
            yield from batches

        for fragment in self._fragments:
            new_fragment = BatchReaderFragment(
                partial(filter_batches, fragment), fragment.schema, fragment._batch_size
            )
            new_fragments.append(new_fragment)

        filtered_dataset = self.__class__(new_fragments)
        filtered_dataset._scan_options = dict(filter=new_filter)

        return filtered_dataset

    def replace_schema(self, schema):
        raise NotImplementedError

    def sort_by(self, sorting, **kwargs):
        raise NotImplementedError

    def join(self, *args, **kwargs):
        raise NotImplementedError

    def join_asof(self, *args, **kwargs):
        raise NotImplementedError

    def to_polars(self):
        """Return a Polars LazyFrame."""
        import polars as pl

        return pl.scan_pyarrow_dataset(self)

    def to_duckdb(self, conn):
        """Return a DuckDB Relation."""
        return conn.from_arrow(self)

    def to_dask(self, find_divisions=False):
        """Return a Dask DataFrame."""
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

    def to_pandas(self):
        """Return a Pandas DataFrame."""
        return self.to_table().to_pandas()

    def to_ipc(self) -> bytes:
        """Serialize as Arrow IPC"""
        s = pa.BufferOutputStream()
        with pa.ipc.new_stream(s, self.schema) as writer:
            writer.write_table(self.to_table())
        buffer = s.getvalue()
        return buffer.to_pybytes()
