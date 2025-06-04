"""
DataSource classes for sequence file formats, including FASTA and FASTQ.
"""

from __future__ import annotations

from typing import Any, Callable, Generator, IO, Self
import pathlib
from typing import IO, Any, Callable, Generator, Self

import pyarrow as pa

from oxbow._core.base import DEFAULT_BATCH_SIZE, DataSource
from oxbow.oxbow import PyFastaScanner, PyFastqScanner


class SequenceFile(DataSource):
    @property
    def _gzi(self) -> str | None:
        return self._gzi_src() if self._gzi_src else None

    def _batchreader_builder(
        self,
    ) -> Callable[[list[str] | None, int], pa.RecordBatchReader]:
        def builder(columns, batch_size):
            scanner = self.scanner()
            scan_kwargs = self._schema_kwargs.copy()
            field_names = scanner.field_names()

            if columns is not None:
                scan_kwargs["fields"] = [col for col in columns if col in field_names]

            if self._regions is not None:
                scan_kwargs["regions"] = self._regions
                scan_kwargs["index"] = self._index
                scan_kwargs["gzi"] = self._gzi
                scan_fn = scanner.scan_query
            else:
                scan_fn = scanner.scan

            stream = scan_fn(**scan_kwargs, batch_size=batch_size)
            return pa.RecordBatchReader.from_stream(
                data=stream,
                schema=pa.schema(stream.schema),
            )

        return builder

    @property
    def _batchreader_builders(
        self,
    ) -> Generator[Callable[[list[str] | None, int], pa.RecordBatchReader]]:
        # Right now, we always produce one fragment, even if multiple regions
        # are requested from a FASTA, as each region corresponds to a single
        # record.
        yield self._batchreader_builder()

    def __init__(
        self,
        source,
        compressed,
        fields,
        regions,
        index,
        gzi,
        batch_size,
    ):
        super().__init__(source, index, batch_size)
        if isinstance(gzi, (str, pathlib.Path)):
            gzi = str(gzi)
            self._gzi_src = lambda: gzi
        elif callable(gzi) or gzi is None:
            self._gzi_src = gzi
        else:
            raise TypeError(
                "`gzi` must be a str, pathlib.Path, or a callable returning "
                "an IO stream"
            )

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._scanner_kwargs = dict(compressed=compressed)
        self._schema_kwargs = dict(fields=fields)

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._src,
            regions=regions,
            index=self._index_src,
            gzi=self._gzi_src,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
            **self._schema_kwargs,
        )


class FastaFile(SequenceFile):
    _scanner_type = PyFastaScanner

    def __init__(
        self,
        source: str | pathlib.Path | Callable[[], IO[Any]],
        compressed: bool = False,
        *,
        fields: list[str] | None = None,
        regions: str | list[str] | None = None,
        index: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
        gzi: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
        batch_size: int = 1,
    ):
        super().__init__(
            source=source,
            compressed=compressed,
            fields=fields,
            regions=regions,
            index=index,
            gzi=gzi,
            batch_size=batch_size,
        )


class FastqFile(SequenceFile):
    _scanner_type = PyFastqScanner

    def __init__(
        self,
        source: str | pathlib.Path | Callable[[], IO[Any]],
        compressed: bool = False,
        *,
        fields: list[str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        super().__init__(
            source=source,
            compressed=compressed,
            fields=fields,
            regions=None,
            index=None,
            gzi=None,
            batch_size=batch_size,
        )

    def regions(self, regions: str | list[str]):
        raise NotImplementedError("FastqFile does not support genomic range queries.")


def from_fasta(
    source: str | pathlib.Path | Callable[[], IO[Any]],
    compression: Literal["infer", "gzip", "bgzf", None] = "infer",
    *,
    fields: list[str] | None = None,
    regions: list[tuple[int, int]] | None = None,
    index: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
    gzi: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
    batch_size: int = 1,
) -> FastaFile:
    """
    Create a FASTA file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the FASTA file, or a callable that opens the file
        as a file-like object.
    compressed : bool, optional
        Whether the source is compressed, by default False.
    fields : list[str], optional
        Names of the fields to project.
    regions : list[tuple[int, int]], optional
        Genomic regions to query.
    index : str, optional
        The FAI index file.
    gzi : str, optional
        The GZI index file for compressed sources.
    batch_size : int, optional [default: 1]
        The size of the batches to read. Since sequences for FASTA files
        can be very long, the default batch size is set to 1 to generate one
        sequence record at a time.

    Returns
    -------
    FastaFile

    See also
    --------
    from_fastq : Create a FASTQ file data source.
    """
    return FastaFile(
        source=source,
        compression=compression,
        fields=fields,
        regions=regions,
        index=index,
        gzi=gzi,
        batch_size=batch_size,
    )


def from_fastq(
    source: str | pathlib.Path | Callable[[], IO[Any]],
    compression: Literal["infer", "gzip", "bgzf", None] = "infer",
    *,
    fields: list[str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> FastqFile:
    """
    Create a FASTQ file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the FASTQ file, or a callable that opens the file
        as a file-like object.
    compressed : bool, optional
        Whether the source is compressed, by default False.
    fields : list[str], optional
        Names of the fields to project.
    batch_size : int, optional
        The size of the batches to read.

    Returns
    -------
    FastqFile

    See also
    --------
    from_fasta : Create a FASTA file data source.
    """
    return FastqFile(
        source=source,
        compression=compression,
        fields=fields,
        batch_size=batch_size,
    )
