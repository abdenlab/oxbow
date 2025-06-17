"""
DataSource classes for sequence file formats, including FASTA and FASTQ.
"""

from __future__ import annotations

import pathlib
from typing import IO, Callable, Generator, Literal

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

import pyarrow as pa

from oxbow._core.base import DEFAULT_BATCH_SIZE, DataSource, prepare_source_and_index
from oxbow.oxbow import PyFastaScanner, PyFastqScanner


class SequenceFile(DataSource):
    @property
    def _gzi(self):
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
        source: str | Callable[[], IO[bytes] | str],
        compressed: bool = False,
        *,
        fields: list[str] | None = None,
        regions: str | list[str] | None = None,
        index: str | Callable[[], IO[bytes] | str] | None = None,
        gzi: str | Callable[[], IO[bytes]] | None = None,
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
        source: str | Callable[[], IO[bytes] | str],
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
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    compression: Literal["infer", "bgzf", "gzip", None] = "infer",
    *,
    fields: list[str] | None = None,
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    gzi: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = 1,
) -> FastaFile:
    """
    Create a FASTA file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the FASTA file, or a callable that opens the file
        as a file-like object.
    compression : Literal["infer", "bgzf", "gzip", None], default: "infer"
        Compression of the source bytestream. If "infer" and ``source`` is a
        URI or path, the file's compression is guessed based on the extension,
        where ".gz" or ".bgz" is interpreted as BGZF. Pass "gzip" to decode
        regular GZIP. If None, the source bytestream is assumed to be
        uncompressed. For more customized decoding, provide a callable
        ``source`` instead.
    fields : list[str], optional
        Specific fields to project. By default, all fields are included.
    regions : list[str], optional
        Provide one or more genomic ranges to slice subsequences as output
        records. Only applicable if an associated index file is available.
    index : str, pathlib.Path, or Callable, optional
        An optional FAI index file associated with the FASTA file. If
        ``source`` is a URI or path and the index file shares the same name
        with a ".fai" extension, the index file is automatically detected.
        If the FASTA file is BGZF-compressed, a GZI index file is also
        required.
    gzi : str, pathlib.Path, or Callable, optional
        An optional GZI index file associated with a BGZF-compressed FASTA
        file. This is required in addition to the FAI index file for random
        access.
    batch_size : int, optional [default: 1]
        The number of records to read in each batch. Since sequences for FASTA
        files can be very long, the default batch size is set to 1 to generate
        one sequence record at a time.

    Returns
    -------
    FastaFile
        A data source object representing the FASTA file.

    See also
    --------
    from_fastq : Create a FASTQ file data source.
    """
    source, index, bgzf_compressed = prepare_source_and_index(
        source, index, compression
    )
    return FastaFile(
        source=source,
        compressed=bgzf_compressed,
        fields=fields,
        regions=regions,
        index=index,
        gzi=gzi,
        batch_size=batch_size,
    )


def from_fastq(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    compression: Literal["infer", "gzip", None] = "infer",
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
    compression : Literal["infer", "gzip", None], default: "infer"
        Compression of the source bytestream. If "infer" and `source` is a URI
        or path, the file's compression is guessed based on the file extension.
        For more custom decoding, provide a callable ``source`` instead.
    fields : list[str], optional
        Specific fields to project. By default, all fields are included.
    batch_size : int, optional [default: 131072]
        The number of records to read in each batch.

    Returns
    -------
    FastqFile
        A data source object representing the FASTQ file.

    Notes
    -----
    Indexed FASTQ files are not supported. Hence, range queries are disallowed
    and files compressed using either regular GZIP or BGZF are decoded using a
    standard GZIP decoder.

    See also
    --------
    from_fasta : Create a FASTA file data source.
    """
    source, _, compressed = prepare_source_and_index(source, None, compression)
    return FastqFile(
        source=source,
        compressed=compressed,
        fields=fields,
        batch_size=batch_size,
    )
