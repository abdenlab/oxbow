"""
DataSource classes for GTF/GFF3 formats.
"""

from __future__ import annotations

import pathlib
from typing import IO, Any, Callable, Generator, Self

import pyarrow as pa

from oxbow._core.base import DEFAULT_BATCH_SIZE, DataSource
from oxbow.oxbow import PyGffScanner, PyGtfScanner


class GxfFile(CompressibleDataSource):
    def __init__(
        self,
        source: str | pathlib.Path | Callable[[], IO[Any]],
        *,
        fields: list[str] | None = None,
        attribute_defs: list[tuple[str, str]] | None = None,
        attribute_scan_rows: int = 1024,
        regions: str | list[str] | None = None,
        index: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
        compression: Literal["infer", "gzip", "bgzf", None] = "infer",
    ):
        super().__init__(source, index, batch_size, compression=compression)

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._scanner_kwargs = dict(compressed=self.compressed)
        if attribute_defs is None:
            attribute_defs = self.scanner().attribute_defs(attribute_scan_rows)
        self._schema_kwargs = dict(fields=fields, attribute_defs=attribute_defs)

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
                if "attributes" not in columns:
                    scan_kwargs["attribute_defs"] = []
                elif scan_kwargs.get("attribute_defs") == []:
                    raise ValueError(
                        "Cannot select `attributes` column if no attribute "
                        "definitions are provided."
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
                yield self._batchreader_builder(
                    scanner.scan_query, scanner.field_names(), region
                )
        else:
            scanner = self.scanner()
            yield self._batchreader_builder(scanner.scan, scanner.field_names())

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._src,
            regions=regions,
            index=self._index_src,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
            **self._schema_kwargs,
        )


class GtfFile(GxfFile):
    _scanner_type = PyGtfScanner


class GffFile(GxfFile):
    _scanner_type = PyGffScanner


def from_gtf(
    source: str | pathlib.Path | Callable[[], IO[Any]],
    compression: Literal["infer", "gzip", "bgzf", None] = "infer",
    *,
    fields: list[str] | None = None,
    attribute_defs: dict | None = None,
    attribute_scan_rows: int = 1024,
    regions: list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> GtfFile:
    """
    Create a GTF file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the GTF file, or a callable that opens the file
        as a file-like object.
    compressed : bool, optional
        Whether the source is compressed, by default False.
    fields : list[str], optional
        Names of the fields to project.
    attribute_defs : dict, optional
        Attribute definitions for the file.
    attribute_scan_rows : int, optional
        Number of rows to scan for attribute definitions, by default 1024.
    regions : list[str], optional
        Genomic regions to query.
    index : str, pathlib.Path, or Callable, optional
        Index file for the GTF file, by default None.
    batch_size : int, optional
        Size of the batch to read.

    Returns
    -------
    GtfFile

    See also
    --------
    from_gff : Create a GFF file data source.
    from_bed : Create a BED file data source.
    from_bigbed : Create a BigBed file data source.
    from_bigwig : Create a BigWig file data source.
    """
    return GtfFile(
        source=source,
        compression=compression,
        fields=fields,
        attribute_defs=attribute_defs,
        attribute_scan_rows=attribute_scan_rows,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )


def from_gff(
    source: str | pathlib.Path | Callable[[], IO[Any]],
    compression: Literal["infer", "gzip", "bgzf", None] = "infer",
    *,
    fields: list[str] | None = None,
    attribute_defs: dict | None = None,
    attribute_scan_rows: int = 1024,
    regions: list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> GffFile:
    """
    Create a GFF3 file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the GFF file, or a callable that opens the file
        as a file-like object.
    compressed : bool, optional
        Whether the source is compressed, by default False.
    fields : list[str], optional
        Names of the fields to project.
    attribute_defs : dict, optional
        Attribute definitions for the file.
    attribute_scan_rows : int, optional
        Number of rows to scan for attribute definitions, by default 1024.
    regions : list[str], optional
        Genomic regions to query.
    index : str, pathlib.Path, or Callable, optional
        Index file for the GFF file, by default None.
    batch_size : int, optional
        Size of the batch to read.

    Returns
    -------
    GffFile

    See also
    --------
    from_gtf : Create a GTF file data source.
    from_bed : Create a BED file data source.
    from_bigbed : Create a BigBed file data source.
    from_bigwig : Create a BigWig file data source.
    """
    return GffFile(
        source=source,
        fields=fields,
        index=index,
        compression=compression,
        regions=regions,
        attribute_defs=attribute_defs,
        attribute_scan_rows=attribute_scan_rows,
        batch_size=batch_size,
    )
