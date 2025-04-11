from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any, Callable, Generator

import pyarrow as pa

from oxbow._core.base import DataFile
from oxbow._filetypes import FileType
from oxbow.oxbow import PyBamScanner, PySamScanner


class AlignmentFile(DataFile):
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
        tag_defs: Any = None,
        tag_scan_rows: int = 1024,
    ) -> None:
        """
        Initialize an AlignmentFile instance.

        Parameters
        ----------
        uri
            The URI or path to the alignment file.
        opener
            A callable to open the file, by default None.
        fields
            Specific fields to extract from the file, by default None.
        index
            Index for the file, by default None.
        compressed
            Whether the file is compressed, by default False.
        regions
            Specific regions to extract from the file, by default None.
        tag_defs
            Definitions for custom tags, by default None.
        tag_scan_rows
            Number of rows to scan for tags, by default 1024.
        """
        super().__init__(uri, opener, fields)
        self._index = index
        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions
        self.__scanner_kwargs = dict(compressed=compressed)
        if not tag_defs:
            tag_defs = self._scanner.tag_defs(tag_scan_rows)
        self.__scan_kwargs = dict(fields=fields, tag_defs=tag_defs)
        self.__schema_kwargs = dict(fields=fields, tag_defs=tag_defs)


class BamFile(AlignmentFile, file_type=FileType.BAM):
    if TYPE_CHECKING:
        _scanner: PyBamScanner

    def __init__(self, *args: Any, compressed: bool = True, **kwargs: Any) -> None:
        """
        Initialize a BamFile instance.

        Parameters
        ----------
        *args
            Positional arguments for the parent class.
        compressed
            Whether the BAM file is compressed, by default True.
        **kwargs
            Keyword arguments for the parent class.
        """
        super().__init__(*args, compressed=compressed, **kwargs)


class SamFile(AlignmentFile, file_type=FileType.SAM):
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
    uri
        The URI or path to the BAM file.
    opener
        A callable to open the file, by default None.
    fields
        Specific fields to extract from the BAM file, by default None.
    index
        Index for the BAM file, by default None.
    compressed
        Whether the BAM file is compressed, by default True.
    regions
        Specific regions to extract from the BAM file, by default None.
    tag_defs
        Definitions for custom tags, by default None.
    tag_scan_rows
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
    uri
        The URI or path to the SAM file.
    opener
        A callable to open the file, by default None.
    fields
        Specific fields to extract from the SAM file, by default None.
    index
        Index for the SAM file, by default None.
    compressed
        Whether the SAM file is compressed, by default False.
    regions
        Specific regions to extract from the SAM file, by default None.
    tag_defs
        Definitions for custom tags, by default None.
    tag_scan_rows
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
