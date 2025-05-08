"""
This module provides classes and utilities for working with PyArrow datasets and fragments.

Classes
-------
BatchReaderFragment
    A fragment that emits RecordBatches from a reproducible source.
BatchReaderDataset
    A PyArrow Dataset composed of one or more BatchReaderFragments.

Constants
---------
DEFAULT_BATCH_SIZE : int
    The default batch size for reading record batches.
DEFAULT_BATCH_READAHEAD : int
    The default number of batches to read ahead.
DEFAULT_FRAGMENT_READAHEAD : int
    The default number of fragments to read ahead.

Type Aliases
------------
RecordBatchIter : Union[pyarrow.RecordBatchReader, Iterator[pyarrow.RecordBatch]]
    A type alias for a PyArrow RecordBatchReader or an iterator of RecordBatches.
"""

from __future__ import annotations

from functools import partial
from itertools import chain
from typing import Callable, Final, Iterator

import pyarrow as pa
import pyarrow.dataset as ds
from pyarrow.dataset import Dataset, Fragment, Scanner

DEFAULT_BATCH_SIZE: Final = 2**17
DEFAULT_BATCH_READAHEAD: Final = 16
DEFAULT_FRAGMENT_READAHEAD: Final = 4

RecordBatchIter = pa.RecordBatchReader | Iterator[pa.RecordBatch]


class BatchReaderFragment(Fragment):
    """
    A Fragment that emits RecordBatches from a reproducible source.

    To provide stateless replay, a new record batch iterator over the same
    records is constructed whenever a scanner is requested.

    Parameters
    ----------
    make_batchreader : Callable[[list[str] | None, int], RecordBatchIter]
        A function that recreates a specific stream of record batches.
    schema : pyarrow.Schema
        The schema of the RecordBatches.
    batch_size : int, optional
        The maximum row count for scanned record batches, by default DEFAULT_BATCH_SIZE.
    partition_expression : pyarrow.dataset.Expression, optional
        A partition expression for the fragment, by default None.

    Attributes
    ----------
    schema : pyarrow.Schema
        The schema of the RecordBatches in the fragment.
    partition_expression : pyarrow.dataset.Expression
        An expression that evaluates to true for all data viewed by this fragment.
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
        make_batchreader : Callable[[list[str] | None, int], RecordBatchIter]
            A function that recreates a specific stream of record batches.
        schema : pyarrow.Schema
            The schema of the RecordBatches.
        batch_size : int, optional
            The maximum row count for scanned record batches, by default DEFAULT_BATCH_SIZE.
        partition_expression : pyarrow.dataset.Expression, optional
            A partition expression for the fragment, by default None.

        Notes
        -----
        The `make_batchreader` function should accept the following arguments
        and return a pyarrow.RecordBatchReader:

        * columns: The columns to project.
        * batch_size: The maximum number of rows per batch.
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
        pyarrow.Schema
            The schema of the RecordBatches in the fragment.
        """
        return self._schema

    @property
    def partition_expression(self) -> ds.Expression:
        """
        Returns
        -------
        pyarrow.dataset.Expression
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
        memory_pool: ds.MemoryPool | None = None,
    ) -> ds.Scanner:
        """
        Build a scan operation against the fragment.

        Parameters
        ----------
        schema : pyarrow.Schema, optional
            The schema to use for scanning. If not specified, uses the Fragment's physical schema.
        columns : list[str], optional
            Names of columns to project. By default, all available columns are projected.
        filter : pyarrow.dataset.Expression, optional
            A filter expression. Scan will return only the rows matching the filter.
        batch_size : int, optional
            The maximum row count for scanned record batches, by default DEFAULT_BATCH_SIZE.
        batch_readahead : int, optional
            The number of batches to read ahead in a file, by default DEFAULT_BATCH_READAHEAD.
        fragment_readahead : int, optional
            The number of fragments/files to read ahead, by default DEFAULT_FRAGMENT_READAHEAD.
        fragment_scan_options : pyarrow.dataset.FragmentScanOptions, optional
            Options specific to a particular scan and fragment type, by default None.
        use_threads : bool, optional
            If enabled, maximum parallelism will be used, by default True.
        memory_pool : pyarrow.dataset.MemoryPool, optional
            For memory allocations, if required. By default, uses the default pool.

        Returns
        -------
        pyarrow.dataset.Scanner
            A scanner object for the fragment.
        """
        # Apply column projection if specified
        schema = schema or self._schema
        batches = self._make_batchreader(columns, batch_size)

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
        batch_size: int = DEFAULT_BATCH_SIZE,
        batch_readahead: int = DEFAULT_BATCH_READAHEAD,
        fragment_readahead: int = DEFAULT_FRAGMENT_READAHEAD,
        fragment_scan_options: ds.FragmentScanOptions | None = None,
        use_threads: bool = True,
        memory_pool: ds.MemoryPool | None = None,
    ) -> Iterator[pa.RecordBatch]:
        """
        Scan and read the fragment as materialized record batches.

        Projections and filters are applied if specified.

        Parameters
        ----------
        schema : pyarrow.Schema, optional
            The schema to use for scanning. If not specified, uses the Fragment's physical schema.
        columns : list[str], optional
            Names of columns to project. By default, all available columns are projected.
        filter : pyarrow.dataset.Expression, optional
            A filter expression. Scan will return only the rows matching the filter.
        batch_size : int, optional
            The maximum row count for scanned record batches, by default DEFAULT_BATCH_SIZE.
        batch_readahead : int, optional
            The number of batches to read ahead in a file, by default DEFAULT_BATCH_READAHEAD.
        fragment_readahead : int, optional
            The number of fragments/files to read ahead, by default DEFAULT_FRAGMENT_READAHEAD.
        fragment_scan_options : pyarrow.dataset.FragmentScanOptions, optional
            Options specific to a particular scan and fragment type, by default None.
        use_threads : bool, optional
            If enabled, maximum parallelism will be used, by default True.
        memory_pool : pyarrow.dataset.MemoryPool, optional
            For memory allocations, if required. By default, uses the default pool.

        Returns
        -------
        Iterator[pyarrow.RecordBatch]
            An iterator of record batches.
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
        """
        Iterate over batches in the fragment.

        Parameters
        ----------
        columns : list[str], optional
            Names of columns to project. By default, all available columns are projected.
        batch_size : int, optional
            The maximum row count for scanned record batches, by default DEFAULT_BATCH_SIZE.

        Returns
        -------
        Iterator[pyarrow.RecordBatch]
            An iterator of record batches.
        """
        return self._make_batchreader(columns, batch_size)

    def __dask_tokenize__(self):
        """
        Tokenize the fragment for Dask.

        Returns
        -------
        tuple
            A tuple representing the tokenized fragment.
        """
        from dask.base import normalize_token

        return (
            normalize_token(self.__class__),
            self._schema,
            self._partition_expression,
            self._batch_size,
        )


class BatchReaderDataset(Dataset):
    """
    A PyArrow Dataset composed of one or more BatchReaderFragments.

    Parameters
    ----------
    fragments : list[BatchReaderFragment]
        The list of fragments.
    partition_expression : pyarrow.dataset.Expression, optional
        A partition expression for the dataset, by default None.

    Attributes
    ----------
    schema : pyarrow.Schema
        The schema of the RecordBatches in the dataset.
    partition_expression : pyarrow.dataset.Expression
        An expression that evaluates to true for all data viewed by this dataset.
    """

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
        partition_expression : pyarrow.dataset.Expression, optional
            A partition expression for the dataset, by default None.
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
        pyarrow.dataset.Expression
            An expression that evaluates to true for all data viewed by this dataset.
        """
        return self._partition_expression

    @property
    def schema(self) -> pa.Schema:
        """
        Returns
        -------
        pyarrow.Schema
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
            An iterator over fragments.

        Notes
        -----
        ``filter`` here is meant to be applied at the fragment level via
        comparison with the partition_expression, not at the row level.
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
            Names of columns to project. By default, all available columns are projected.
        filter : pyarrow.dataset.Expression, optional
            A filter expression. Scan will return only the rows matching the filter.
        batch_size : int, optional
            The maximum row count for scanned record batches, by default DEFAULT_BATCH_SIZE.
        batch_readahead : int, optional
            The number of batches to read ahead in a file, by default DEFAULT_BATCH_READAHEAD.
        fragment_readahead : int, optional
            The number of fragments/files to read ahead, by default DEFAULT_FRAGMENT_READAHEAD.
        fragment_scan_options : pyarrow.dataset.FragmentScanOptions, optional
            Options specific to a particular scan and fragment type, by default None.
        use_threads : bool, optional
            If enabled, maximum parallelism will be used, by default True.
        memory_pool : pyarrow.dataset.MemoryPool, optional
            For memory allocations, if required. By default, uses the default pool.

        Returns
        -------
        pyarrow.dataset.Scanner
            A scanner object for the dataset.
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
            Names of columns to project. By default, all available columns are projected.
        filter : pyarrow.dataset.Expression, optional
            A filter expression. Scan will return only the rows matching the filter.
        batch_size : int, optional
            The maximum row count for scanned record batches, by default DEFAULT_BATCH_SIZE.
        batch_readahead : int, optional
            The number of batches to read ahead in a file, by default DEFAULT_BATCH_READAHEAD.
        fragment_readahead : int, optional
            The number of fragments/files to read ahead, by default DEFAULT_FRAGMENT_READAHEAD.
        fragment_scan_options : pyarrow.dataset.FragmentScanOptions, optional
            Options specific to a particular scan and fragment type, by default None.
        use_threads : bool, optional
            If enabled, maximum parallelism will be used, by default True.
        memory_pool : pyarrow.dataset.MemoryPool, optional
            For memory allocations, if required. By default, uses the default pool.

        Returns
        -------
        Iterator[pyarrow.RecordBatch]
            An iterator of record batches.
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
        """
        Iterate over batches in the dataset.

        Parameters
        ----------
        columns : list[str], optional
            Names of columns to project. By default, all available columns are projected.
        batch_size : int, optional
            The maximum row count for scanned record batches, by default DEFAULT_BATCH_SIZE.

        Returns
        -------
        Iterator[pyarrow.RecordBatch]
            An iterator of record batches.
        """
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
        expression : pyarrow.dataset.Expression
            The filter expression.

        Returns
        -------
        BatchReaderDataset
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
        """
        Replace the schema of the dataset.

        Parameters
        ----------
        schema : pyarrow.Schema
            The new schema.

        Raises
        ------
        NotImplementedError
            This method is not yet implemented.
        """
        raise NotImplementedError

    def sort_by(self, sorting, **kwargs):
        """
        Sort the dataset by the specified columns.

        Parameters
        ----------
        sorting : list[str]
            The columns to sort by.

        Raises
        ------
        NotImplementedError
            This method is not yet implemented.
        """
        raise NotImplementedError

    def join(self, *args, **kwargs):
        """
        Perform a join operation on the dataset.

        Raises
        ------
        NotImplementedError
            This method is not yet implemented.
        """
        raise NotImplementedError

    def join_asof(self, *args, **kwargs):
        """
        Perform an as-of join operation on the dataset.

        Raises
        ------
        NotImplementedError
            This method is not yet implemented.
        """
        raise NotImplementedError

    def to_polars(self):
        """
        Convert the dataset to a Polars LazyFrame.

        Returns
        -------
        polars.LazyFrame
            A Polars LazyFrame representation of the dataset.
        """
        import polars as pl

        return pl.scan_pyarrow_dataset(self)

    def to_duckdb(self, conn):
        """
        Convert the dataset to a DuckDB Relation.

        Parameters
        ----------
        conn : duckdb.DuckDBPyConnection
            The DuckDB connection.

        Returns
        -------
        duckdb.DuckDBPyRelation
            A DuckDB Relation representation of the dataset.
        """
        return conn.from_arrow(self)

    def to_dask(self, find_divisions=False):
        """
        Convert the dataset to a Dask DataFrame.

        Parameters
        ----------
        find_divisions : bool, optional
            If True, find divisions for the Dask DataFrame, by default False.

        Returns
        -------
        dask.dataframe.DataFrame
            A Dask DataFrame representation of the dataset.
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

    def to_pandas(self):
        """
        Convert the dataset to a Pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
            A Pandas DataFrame representation of the dataset.
        """
        return self.to_table().to_pandas()

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
            writer.write_table(self.to_table())
        buffer = s.getvalue()
        return buffer.to_pybytes()
