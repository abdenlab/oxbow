"""
This module defines classes and functions for working with alignment files, including BAM and SAM formats.

Classes
-------
AlignmentFile
    Base class for alignment files.
BamFile
    Class for handling BAM files.
SamFile
    Class for handling SAM files.

Functions
---------
from_bam(uri, opener, fields, index, compressed, regions, tag_defs, tag_scan_rows)
    Create a BamFile instance from a BAM file.
from_sam(uri, opener, fields, index, compressed, regions, tag_defs, tag_scan_rows)
    Create a SamFile instance from a SAM file.
"""

from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any, Callable, Generator

import pyarrow as pa

from oxbow._core.base import DataFile
from oxbow._filetypes import FileType
from oxbow.oxbow import PyBamScanner, PySamScanner


class AlignmentFile(DataFile):
    """
    Base class for alignment files.

    This class provides common functionality for handling BAM and SAM files.

    Functions
    ---------
    from_bam(uri, opener, fields, index, compressed, regions, tag_defs, tag_scan_rows)
        Create a BamFile instance from a BAM file.
    from_sam(uri, opener, fields, index, compressed, regions, tag_defs, tag_scan_rows)
        Create a SamFile instance from a SAM file.
    """

    if TYPE_CHECKING:
        _scanner: PyBamScanner | PySamScanner
        from_bam: BamFile
        from_sam: SamFile

    @property
    def _scan_kwargs(self) -> dict[str, Any]:
        return self.__scan_kwargs

    @property
    def _scanner_kwargs(self) -> dict[str, Any]:
        return self.__scanner_kwargs

    @property
    def _schema_kwargs(self) -> dict[str, Any]:
        return self.__schema_kwargs

    @property
    def _batch_readers(
        self,
    ) -> Generator[Callable[[list[str] | None, int], pa.RecordBatchReader]]:
        if self._regions:
            for region in self._regions:
                if region == "*":
                    stream = partial(
                        self._scanner.scan_unmapped,
                        index=self._index,
                        **self._scan_kwargs,
                    )
                else:
                    stream = partial(
                        self._scanner.scan_query,
                        region=region,
                        index=self._index,
                        **self._scan_kwargs,
                    )
                yield lambda fields, batch_size: pa.RecordBatchReader.from_stream(
                    data=stream(fields=fields, batch_size=batch_size),
                    schema=self.schema,
                )
        else:
            stream = partial(
                self._scanner.scan,
                **self._scan_kwargs,
            )
            yield lambda fields, batch_size: pa.RecordBatchReader.from_stream(
                data=stream(fields=fields, batch_size=batch_size),
                schema=self.schema,
            )

    def __init__(
        self,
        uri: str,
        opener: Callable | None = None,
        fields: list[str] | None = None,
        index: Any = None,
        compressed: bool = False,
        regions: str | list[str] | None = None,
        tag_defs: list[tuple[str, str]] | None = None,
        tag_scan_rows: int = 1024,
    ) -> None:
        """
        Initialize an AlignmentFile instance.

        Parameters
        ----------
        uri : str
            The URI or path to the alignment file.
        opener : Callable, optional
            A callable to open the file, by default None.
        fields : list[str], optional
            Specific fields to extract from the file, by default None.
        index : Any, optional
            Index for the file, by default None.
        compressed : bool, optional
            Whether the file is compressed, by default False.
        regions : str | list[str], optional
            Specific regions to extract from the file, by default None.
        tag_defs : Any, optional
            Definitions for custom tags. If None, the scanner will scan the
            first `tag_scan_rows` rows to find tag definitions.
        tag_scan_rows : int, optional
            Number of rows to scan for tag definitions, if tag_defs is None.
            By default 1024.
        """
        super().__init__(uri, opener, fields)
        self._index = index
        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions
        self.__scanner_kwargs = dict(compressed=compressed)
        if tag_defs is None:
            tag_defs = self._scanner.tag_defs(tag_scan_rows)
        self.__scan_kwargs = dict(fields=fields, tag_defs=tag_defs)
        self.__schema_kwargs = dict(fields=fields, tag_defs=tag_defs)


class BamFile(AlignmentFile, file_type=FileType.BAM):
    """
    Class for handling BAM files.

    Parameters
    ----------
    *args : Any
        Positional arguments for the parent class.
    compressed : bool, optional
        Whether the BAM file is compressed, by default True.
    **kwargs : Any
        Keyword arguments for the parent class.
    """

    if TYPE_CHECKING:
        _scanner: PyBamScanner

    def __init__(self, *args: Any, compressed: bool = True, **kwargs: Any) -> None:
        """
        Initialize a BamFile instance.

        Parameters
        ----------
        *args : Any
            Positional arguments for the parent class.
        compressed : bool, optional
            Whether the BAM file is compressed, by default True.
        **kwargs : Any
            Keyword arguments for the parent class.
        """
        super().__init__(*args, compressed=compressed, **kwargs)


class SamFile(AlignmentFile, file_type=FileType.SAM):
    """
    Class for handling SAM files.
    """

    if TYPE_CHECKING:
        _scanner: PySamScanner


def from_bam(
    uri: str,
    opener: Callable | None = None,
    fields: list[str] | None = None,
    index: str | None = None,
    compressed: bool = True,
    regions: str | list[str] | None = None,
    tag_defs: Any = None,
    tag_scan_rows: int = 1024,
) -> BamFile:
    """
    Create a BamFile instance from a BAM file.

    Parameters
    ----------
    uri : str
        The URI or path to the BAM file.
    opener : Callable, optional
        A callable to open the file, by default None.
    fields : list[str], optional
        Specific fields to extract from the BAM file, by default None.
    index : str, optional
        Index for the BAM file, by default None.
    compressed : bool, optional
        Whether the BAM file is compressed, by default True.
    regions : str | list[str], optional
        Specific regions to extract from the BAM file, by default None.
    tag_defs : Any, optional
        Definitions for custom tags, by default None.
    tag_scan_rows : int, optional
        Number of rows to scan for tags, by default 1024.

    Returns
    -------
    BamFile
        An instance of BamFile representing the parsed BAM file.

    Notes
    -----
    BAM (Binary Alignment Map) is a compressed binary representation of SAM files.
    """
    return BamFile(
        uri=uri,
        opener=opener,
        fields=fields,
        index=index,
        compressed=compressed,
        regions=regions,
        tag_defs=tag_defs,
        tag_scan_rows=tag_scan_rows,
    )


def from_sam(
    uri: str,
    opener: Callable | None = None,
    fields: list[str] | None = None,
    index: Any = None,
    compressed: bool = False,
    regions: str | list[str] | None = None,
    tag_defs: Any = None,
    tag_scan_rows: int = 1024,
) -> SamFile:
    """
    Create a SamFile instance from a SAM file.

    Parameters
    ----------
    uri : str
        The URI or path to the SAM file.
    opener : Callable, optional
        A callable to open the file, by default None.
    fields : list[str], optional
        Specific fields to extract from the SAM file, by default None.
    index : Any, optional
        Index for the SAM file, by default None.
    compressed : bool, optional
        Whether the SAM file is compressed, by default False.
    regions : str | list[str], optional
        Specific regions to extract from the SAM file, by default None.
    tag_defs : Any, optional
        Definitions for custom tags, by default None.
    tag_scan_rows : int, optional
        Number of rows to scan for tags, by default 1024.

    Returns
    -------
    SamFile
        An instance of SamFile representing the parsed SAM file.

    Notes
    -----
    SAM (Sequence Alignment Map) is a widely used text-based format for storing
    biological sequences aligned to a reference sequence.
    """
    return SamFile(
        uri=uri,
        opener=opener,
        fields=fields,
        index=index,
        compressed=compressed,
        regions=regions,
        tag_defs=tag_defs,
        tag_scan_rows=tag_scan_rows,
    )
