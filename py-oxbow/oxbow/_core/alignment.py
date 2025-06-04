"""
DataSource classes for htslib alignment formats.
"""

from __future__ import annotations

import pathlib
from typing import IO, Any, Callable, Generator, Self

import pyarrow as pa

from oxbow._core.base import DEFAULT_BATCH_SIZE, DataSource
from oxbow.oxbow import PyBamScanner, PySamScanner


class AlignmentFile(CompressibleDataSource):
    def _batchreader_builder(
        self,
        scan_fn: Callable,
        field_names: list[str],
        region: str | None = None,
    ) -> Callable[[list[str] | None, int], pa.RecordBatchReader]:
        def builder(columns, batch_size):
            scan_kwargs = self._schema_kwargs.copy()

            if columns is not None:
                scan_kwargs["fields"] = [col for col in columns if col in field_names]
                if "tags" not in columns:
                    scan_kwargs["tag_defs"] = []
                elif scan_kwargs.get("tag_defs") == []:
                    raise ValueError(
                        "Cannot select `tags` column if no tag definitions are "
                        "provided."
                    )

            if region is not None:
                scan_kwargs["region"] = region
                scan_kwargs["index"] = self._index

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
        if self._regions:
            for region in self._regions:
                scanner = self.scanner()
                if region == "*":
                    yield self._batchreader_builder(
                        scanner.scan_unmapped, scanner.field_names()
                    )
                else:
                    yield self._batchreader_builder(
                        scanner.scan_query, scanner.field_names(), region
                    )
        else:
            scanner = self.scanner()
            yield self._batchreader_builder(scanner.scan, scanner.field_names())

    def __init__(
        self,
        source: str | pathlib.Path | Callable[[], IO[Any]],
        compression: Literal["infer", "gzip", "bgzf", None] = "infer",
        *,
        fields: list[str] | None = None,
        tag_defs: list[tuple[str, str]] | None = None,
        tag_scan_rows: int = 1024,
        regions: str | list[str] | None = None,
        index: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        super().__init__(source, index, batch_size, compression=compression)

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._scanner_kwargs = dict(compressed=self.compressed)
        if tag_defs is None:
            tag_defs = self.scanner().tag_defs(tag_scan_rows)
        self._schema_kwargs = dict(fields=fields, tag_defs=tag_defs)

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._src,
            regions=regions,
            index=self._index_src,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
            **self._schema_kwargs,
        )


class SamFile(AlignmentFile):
    _scanner_type = PySamScanner


class BamFile(AlignmentFile):
    _scanner_type = PyBamScanner

    def __init__(self, *args, compression: Literal["bgzf", None] = "bgzf", **kwargs):
        super().__init__(*args, compression=compression, **kwargs)


def from_sam(
    source: str | pathlib.Path | Callable[[], IO[Any]],
    compressed: bool = False,
    *,
    fields: list[str] | None = None,
    tag_defs: Any = None,
    tag_scan_rows: int = 1024,
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> SamFile:
    """
    Create a SAM (Sequence Alignment Map) file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the SAM file, or a callable that opens the file
        as a file-like object.
    compressed : bool, optional
        Whether the SAM file is compressed, by default False.
    fields : list[str], optional
        Specific fields to extract from the SAM file, by default None.
    tag_defs : Any, optional
        Definitions for custom tags, by default None.
    tag_scan_rows : int, optional
        Number of rows to scan for tags, by default 1024.
    regions : str | list[str], optional
        Specific regions to extract from the SAM file, by default None.
    index : str, pathlib.Path, or Callable, optional
        Index file for the SAM file, by default None.
    batch_size : int, optional
        The size of the batches to read.

    Returns
    -------
    SamFile
        A data source object representing the SAM file.

    Notes
    -----
    SAM is a widely used text-based format for storing biological sequences
    aligned to a reference sequence.

    See also
    --------
    from_bam : Create a BAM file data source.
    """
    return SamFile(
        source=source,
        compressed=compressed,
        fields=fields,
        tag_defs=tag_defs,
        tag_scan_rows=tag_scan_rows,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )


def from_bam(
    source: str | pathlib.Path | Callable[[], IO[Any]],
    compressed: bool = True,
    *,
    fields: list[str] | None = None,
    tag_defs: Any = None,
    tag_scan_rows: int = 1024,
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> BamFile:
    """
    Create a BAM (Binary Alignment Map) file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the BAM file, or a callable that opens the file
        as a file-like object.
    compressed : bool, optional
        Whether the BAM file is compressed, by default True.
    fields : list[str], optional
        Specific fields to extract from the BAM file, by default None.
    tag_defs : Any, optional
        Definitions for custom tags, by default None.
    tag_scan_rows : int, optional
        Number of rows to scan for tags, by default 1024.
    regions : str | list[str], optional
        Specific regions to extract from the BAM file, by default None.
    index : str, pathlib.Path, or Callable, optional
        Index file for the BAM file, by default None.
    batch_size : int, optional
        The size of the batches to read.

    Returns
    -------
    BamFile
         A data source object representing the BAM file.

    Notes
    -----
    BAM is a compressed binary representation of SAM files.

    See also
    --------
    from_sam : Create a SAM file data source.
    """
    return BamFile(
        source=source,
        compressed=compressed,
        fields=fields,
        tag_defs=tag_defs,
        tag_scan_rows=tag_scan_rows,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )
