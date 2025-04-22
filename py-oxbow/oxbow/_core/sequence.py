"""
This module defines classes and functions for working with sequence files, including FASTA and FASTQ formats.

Classes
-------
SequenceFile
    Base class for sequence files.
FastaFile
    Class for handling FASTA files.
FastqFile
    Class for handling FASTQ files.

Functions
---------
from_fasta(uri, opener, fields, index, gzi, compressed, regions)
    Create a FastaFile instance from a FASTA file.
from_fastq(uri, opener, fields, index, gzi, compressed)
    Create a FastqFile instance from a FASTQ file.
"""

from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any, Callable, Generator

import pyarrow as pa

from oxbow._core.base import DataFile
from oxbow._filetypes import FileType
from oxbow._pyarrow import BatchReaderDataset
from oxbow.oxbow import PyFastaScanner, PyFastqScanner


class SequenceFile(DataFile):
    """
    Base class for sequence files.

    This class provides common functionality for handling FASTA and FASTQ files.

    Functions
    ---------
    from_fasta(uri, opener, fields, index, gzi, compressed, regions)
        Create a FastaFile instance from a FASTA file.
    from_fastq(uri, opener, fields, index, gzi, compressed)
        Create a FastqFile instance from a FASTQ file.
    """
    if TYPE_CHECKING:
        from_fasta: FastaFile
        from_fastq: FastqFile
        _scanner: PyFastaScanner | PyFastqScanner

    @property
    def _scan_kwargs(self) -> dict[str, Any]:
        return self.__scan_kwargs

    @property
    def _scanner_kwargs(self) -> dict[str, Any]:
        return self.__scanner_kwargs

    @property
    def _schema_kwargs(self) -> dict[str, Any]:
        return self.__schema_kwargs

    def __init__(
        self,
        uri,
        opener=None,
        fields=None,
        index=None,
        gzi=None,
        compressed=False,
    ):
        self._index = index
        self._gzi = gzi
        self.__schema_kwargs = dict(fields=fields)
        self.__scan_kwargs = dict(fields=fields)
        self.__scanner_kwargs = dict(compressed=compressed)
        super().__init__(uri, opener, fields)

    def select(self, *args, **kwargs) -> BatchReaderDataset:
        return super().select(*args, gzi=self._gzi, index=self._index, **kwargs)


class FastaFile(SequenceFile, file_type=FileType.FASTA):
    """
    Class for handling FASTA files.

    Parameters
    ----------
    uri : str
        The URI of the FASTA file.
    opener : Callable, optional
        A callable to open the file.
    fields : list[str], optional
        Names of the fields to project.
    index : str, optional
        Path to the FAI index file.
    gzi : str, optional
        Path to the GZI index file for compressed sources.
    compressed : bool, optional
        Whether the source is compressed, by default False.
    regions : list[tuple[int, int]], optional
        Genomic regions to query.
    """
    if TYPE_CHECKING:
        _scanner: FileType.FASTA.value

    @property
    def _batch_readers(
        self,
    ) -> Generator[Callable[[list[str] | None, int], pa.RecordBatchReader]]:
        if self._regions:
            stream = partial(
                self._scanner.scan_query,
                regions=self._regions,
                index=self._index,
                gzi=self._gzi,
            )
            yield lambda fields, batch_size: pa.RecordBatchReader.from_stream(
                data=stream(fields=fields, batch_size=batch_size),
                schema=self.schema,
            )
        else:
            stream = self._scanner.scan
            yield lambda fields, batch_size: pa.RecordBatchReader.from_stream(
                data=stream(fields=fields, batch_size=batch_size),
                schema=self.schema,
            )

    def __init__(
        self,
        uri,
        opener=None,
        fields=None,
        index=None,
        gzi=None,
        compressed=False,
        regions=None,
    ):
        self._regions = regions
        super().__init__(uri, opener, fields, index, gzi, compressed)

    def select(self, *args, **kwargs) -> BatchReaderDataset:
        return super().select(*args, regions=self._regions, **kwargs)


class FastqFile(SequenceFile, file_type=FileType.FASTQ):
    """
    Class for handling FASTQ files.

    Parameters
    ----------
    uri : str
        The URI of the FASTQ file.
    opener : Callable, optional
        A callable to open the file.
    fields : list[str], optional
        Names of the fields to project.
    index : str, optional
        Path to the index file.
    gzi : str, optional
        Path to the GZI index file for compressed sources.
    compressed : bool, optional
        Whether the source is compressed, by default False.
    """
    if TYPE_CHECKING:
        _scanner: FileType.FASTQ.value

    @property
    def _batch_readers(
        self,
    ) -> Generator[Callable[[list[str] | None, int], pa.RecordBatchReader]]:
        stream = self._scanner.scan
        yield lambda fields, batch_size: pa.RecordBatchReader.from_stream(
            data=stream(fields=fields, batch_size=batch_size),
            schema=self.schema,
        )


def from_fasta(
    uri: str,
    opener: Callable | None = None,
    fields: list[str] | None = None,
    index: str | None = None,
    gzi: str | None = None,
    compressed: bool = False,
    regions: list[tuple[int, int]] | None = None,
) -> FastaFile:
    """
    Create a FastaFile instance from a FASTA file.

    Parameters
    ----------
    uri : str
        The URI of the FASTA file.
    opener : Callable, optional
        A callable to open the file.
    fields : list[str], optional
        Names of the fields to project.
    index : str, optional
        Path to the FAI index file.
    gzi : str, optional
        Path to the GZI index file for compressed sources.
    compressed : bool, optional
        Whether the source is compressed, by default False.
    regions : list[tuple[int, int]], optional
        Genomic regions to query.

    Returns
    -------
    FastaFile
        An instance of the FastaFile class.
    """
    return FastaFile(
        uri=uri,
        opener=opener,
        fields=fields,
        index=index,
        gzi=gzi,
        compressed=compressed,
        regions=regions,
    )


def from_fastq(
    uri: str,
    opener: Callable | None = None,
    fields: list[str] | None = None,
    index: str | None = None,
    gzi: str | None = None,
    compressed: bool = False,
) -> FastqFile:
    """
    Create a FastqFile instance from a FASTQ file.

    Parameters
    ----------
    uri : str
        The URI of the FASTQ file.
    opener : Callable, optional
        A callable to open the file.
    fields : list[str], optional
        Names of the fields to project.
    index : str, optional
        Path to the index file.
    gzi : str, optional
        Path to the GZI index file for compressed sources.
    compressed : bool, optional
        Whether the source is compressed, by default False.

    Returns
    -------
    FastqFile
        An instance of the FastqFile class.
    """
    return FastqFile(
        uri=uri,
        opener=opener,
        fields=fields,
        index=index,
        gzi=gzi,
        compressed=compressed,
    )
