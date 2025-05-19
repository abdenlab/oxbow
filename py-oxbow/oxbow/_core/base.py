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
    _registry = {}

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

    def _make_batch_reader(self, stream):
        def _batch_reader(fields, batch_size):
            return pa.RecordBatchReader.from_stream(
                data=stream(fields=fields, batch_size=batch_size),
                schema=self.schema,
            )

        return _batch_reader

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
        *,
        batch_size: int = DEFAULT_BATCH_SIZE,
        **kwargs,
    ) -> BatchReaderDataset:
        """
        Select a subset of the data file.

        This method creates a new instance of the data file with the same
        parameters, applies any overrides specified by keyword arguments,
        and returns it as a dataset.

        Parameters
        ----------
        batch_size
            The size of each batch to be read. Default is DEFAULT_BATCH_SIZE.
        **kwargs
            Additional keyword arguments to override the current instance's
            parameters.

        Returns
        -------
        BatchReaderDataset
            A dataset containing the selected subset of the data file.
        """
        return type(self)(
            self._uri,
            self._opener,
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
        Generate record batches from the data file.

        This method iterates over the data file and yields record batches
        of the specified size.

        Parameters
        ----------
        batch_size
            The size of each batch to be read. Default is DEFAULT_BATCH_SIZE.

        Yields
        ------
        pa.RecordBatch
            A record batch from the data file.
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
        Get fragments of the data file.

        Fragments represent parts of the data file that can be processed
        independently. Each fragment is associated with a schema and batch size.

        Parameters
        ----------
        batch_size
            The size of each batch to be read. Default is DEFAULT_BATCH_SIZE.

        Returns
        -------
        list of BatchReaderFragment
            A list of fragments representing parts of the data file.
        """
        return [
            BatchReaderFragment(r, self.schema, batch_size=batch_size)
            for r in self._batch_readers
        ]

    def dataset(self, *, batch_size=DEFAULT_BATCH_SIZE) -> BatchReaderDataset:
        """
        Convert the data file into a dataset.

        A dataset is a collection of fragments that can be processed
        as a single logical entity.

        Parameters
        ----------
        batch_size
            The size of each batch to be read. Default is DEFAULT_BATCH_SIZE.

        Returns
        -------
        BatchReaderDataset
            A dataset representation of the data file.
        """
        return BatchReaderDataset(self.fragments(batch_size=batch_size))
