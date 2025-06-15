import io
import pathlib
import pytest
import itertools
from unittest.mock import MagicMock

import pyarrow as pa
import oxbow.core as ox
import polars as pl

from oxbow._core.base import prepare_source_and_index

from oxbow._pyarrow import (
    BatchReaderDataset,
    BatchReaderFragment,
)


class TestDataSourceInitialization:
    source = "test.src"
    index = "test.idx"

    def test_init_string_source(self):
        source = "test.src"
        ds = ox.DataSource(source)
        assert ds._src == source
        assert ds._source == source

    def test_init_callable_source(self):
        ds = ox.DataSource(lambda: self.source)
        assert callable(ds._src)
        assert ds._source == self.source

    def test_init_with_index(self):
        ds = ox.DataSource(self.source, index=self.index)
        assert ds._index_src == self.index
        assert ds._index == self.index

    def test_init_with_callable_index(self):
        ds = ox.DataSource(self.source, index=lambda: self.index)
        assert callable(ds._index_src)
        assert ds._index == self.index

    def test_init_with_batch_size(self):
        ds = ox.DataSource(self.source, batch_size=42)
        assert ds._batch_size == 42

    def test_init_no_args(self):
        with pytest.raises(TypeError):
            ox.DataSource()


def test_scanner():
    mock_scanner_type = MagicMock()
    source = "test.src"
    scanner_kwargs = {"foo": "bar", "baz": 2}

    ds = ox.DataSource(source)
    ds._scanner_type = mock_scanner_type
    ds._scanner_kwargs = scanner_kwargs

    scanner = ds.scanner()
    mock_scanner_type.assert_called_once_with(source, **scanner_kwargs)
    assert scanner == mock_scanner_type.return_value


def test_schema_and_columns():
    source = "test.src"
    expected_schema = pa.schema([("name", pa.string())])
    expected_schema_kwargs = {"foo": "bar", "baz": 2}

    ds = ox.DataSource(source)
    ds._schema_kwargs = expected_schema_kwargs

    mocked_scanner_instance = MagicMock()
    mocked_scanner_instance.schema.return_value = expected_schema
    ds.scanner = lambda: mocked_scanner_instance

    assert ds.schema == expected_schema
    mocked_scanner_instance.schema.assert_called_once_with(**ds._schema_kwargs)

    assert ds.columns == ["name"]


def test_batches():
    batch_size = 10
    columns = ["name", "age"]

    class MockRecordBatchReader:
        def __init__(self, n_batches):
            self.count = 0
            self.n_batches = n_batches

        def read_next_batch(self):
            self.count += 1
            if self.count <= self.n_batches:
                return MagicMock()
            else:
                raise StopIteration

    mock_builder_one = MagicMock()
    mock_builder_one.return_value = MockRecordBatchReader(n_batches=1)

    mock_builder_two = MagicMock()
    mock_builder_two.return_value = MockRecordBatchReader(n_batches=3)

    class MockDataSource(ox.DataSource):
        @property
        def _batchreader_builders(self):
            return [mock_builder_one, mock_builder_two]

        @property
        def columns(self):
            return columns

        def scanner(self):
            pass

    ds = MockDataSource("test.src", batch_size=batch_size)

    assert len(list(ds.batches())) == 4
    mock_builder_one.assert_called_with(columns, batch_size)
    mock_builder_two.assert_called_with(columns, batch_size)


def test_fragments_and_dataset():
    batch_size = 9

    mock_builder_one = MagicMock()
    mock_builder_two = MagicMock()
    mock_builder_three = MagicMock()

    class MockDataSource(ox.DataSource):
        @property
        def _batchreader_builders(self):
            return [mock_builder_one, mock_builder_two, mock_builder_three]

        @property
        def schema(self):
            return MagicMock()

    ds = MockDataSource("test.src", batch_size=batch_size)

    assert len(ds.fragments()) == 3
    assert isinstance(ds.dataset(), BatchReaderDataset)


def test_to_pandas():
    import pandas as pd

    class MockDataSource(ox.DataSource):
        def dataset(self):
            return BatchReaderDataset(
                [
                    BatchReaderFragment(
                        MagicMock(), pa.schema([("name", pa.string())]), batch_size=10
                    )
                ]
            )

    ds = MockDataSource("test.src", batch_size=10)
    assert isinstance(ds.to_pandas(), pd.DataFrame)


def test_to_polars_eager():
    class MockDataSource(ox.DataSource):
        def batches(self):
            return [pa.RecordBatch.from_pydict({"name": ["test"]})]

    ds = MockDataSource("test.src")
    assert isinstance(ds.to_polars(lazy=False), pl.DataFrame)


def test_to_polars_lazy():
    class MockDataSource(ox.DataSource):
        def dataset(self):
            return BatchReaderDataset(
                [
                    BatchReaderFragment(
                        MagicMock(), pa.schema([("name", pa.string())]), batch_size=10
                    )
                ]
            )

    ds = MockDataSource("test.src", batch_size=10)
    assert isinstance(ds.to_polars(lazy=True), pl.LazyFrame)


def test_to_dask():
    import dask.dataframe as dd

    class MockDataSource(ox.DataSource):
        @property
        def schema(self):
            return pa.schema([("name", pa.string())])

        def fragments(self):
            return [
                BatchReaderFragment(
                    MagicMock(), pa.schema([("name", pa.string())]), batch_size=10
                )
            ]

    ds = MockDataSource("test.src", batch_size=10)
    assert isinstance(ds.to_dask(), dd.DataFrame)


def test_to_duckdb():
    import duckdb

    class MockDataSource(ox.DataSource):
        def dataset(self):
            return BatchReaderDataset(
                [
                    BatchReaderFragment(
                        MagicMock(), pa.schema([("name", pa.string())]), batch_size=10
                    )
                ]
            )

    conn = duckdb.connect()
    ds = MockDataSource("test.src", batch_size=10)
    assert isinstance(ds.to_duckdb(conn), duckdb.DuckDBPyRelation)


def test_to_ipc():
    import pyarrow.ipc as ipc

    class MockDataSource(ox.DataSource):
        @property
        def schema(self):
            return pa.schema([("name", pa.string())])

        def dataset(self):
            return BatchReaderDataset(
                [
                    BatchReaderFragment(
                        MagicMock(), pa.schema([("name", pa.string())]), batch_size=10
                    )
                ]
            )

    ds = MockDataSource("test.src", batch_size=10)
    assert isinstance(ds.to_ipc(), bytes)
    assert ipc.open_stream(ds.to_ipc()).schema == pa.schema([("name", pa.string())])


def dummy_bytes_io():
    return io.BytesIO(b"test")


@pytest.mark.parametrize(
    "source_arg, index_arg, compression_arg",
    itertools.product(
        [
            "file.src",
            pathlib.Path("file.src"),
            "file.src.gz",
            "file.src.bgz",
            dummy_bytes_io,
        ],
        ["file.idx", pathlib.Path("file.idx"), dummy_bytes_io, None],
        ["infer", "bgzf", "gzip", None],
    ),
)
def test_prepare_source_and_index_exhaustive(source_arg, index_arg, compression_arg):
    # Special case for callable sources and gzip compression
    if callable(source_arg) and compression_arg == "gzip":
        with pytest.raises(ValueError):
            prepare_source_and_index(source_arg, index_arg, compression_arg)
        return

    actual_src, actual_index_src, actual_bgzf_compressed = prepare_source_and_index(
        source_arg, index_arg, compression_arg
    )

    assert isinstance(actual_src, (str, pathlib.Path, io.BytesIO, type(lambda: None)))
    assert isinstance(
        actual_index_src,
        (str, pathlib.Path, io.BytesIO, type(lambda: None), type(None)),
    )
    assert isinstance(actual_bgzf_compressed, bool)


@pytest.mark.parametrize(
    "url",
    [
        "https://example.com/file.src",
        "http://example.com/file.src",
        "s3://bucket/file.src",
        "ftp://server/file.src",
        "file://file.src",
    ],
)
def test_prepare_source_and_index_with_fsspec(url):
    actual_src, actual_index_src, _ = prepare_source_and_index(url, url, None)

    assert isinstance(actual_src, type(lambda: None))
    assert isinstance(actual_index_src, type(lambda: None))


@pytest.mark.parametrize("compressed_file", ["file.src.gz", "file.src.bgz"])
def test_prepare_source_and_index_infer_compression(compressed_file):
    _, _, actual_bgzf_compressed = prepare_source_and_index(
        compressed_file, compressed_file, "infer"
    )

    assert actual_bgzf_compressed is True


def test_prepare_source_and_index_gzip():
    actual_src, _, actual_bgzf_compressed = prepare_source_and_index(
        "file.src", None, "gzip"
    )
    assert isinstance(actual_src, type(lambda: None))
    assert actual_bgzf_compressed is False
