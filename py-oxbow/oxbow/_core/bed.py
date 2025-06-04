"""
DataSource classes for the BED family of formats.
"""

from __future__ import annotations

import pathlib
from typing import IO, Any, Callable, Generator, Self

import pyarrow as pa

from oxbow._core.base import DEFAULT_BATCH_SIZE, DataSource
from oxbow.oxbow import PyBedScanner


class BedFile(DataSource):
    _scanner_type = PyBedScanner

    def _batchreader_builder(
        self,
        scan_fn: Callable,
        field_names: list[str],
        region: str | None = None,
    ) -> Callable[[list[str] | None, int], pa.RecordBatchReader]:
        def builder(columns, batch_size):
            scan_kwargs = self._schema_kwargs.copy()

            # TODO: Subsetting BED fields is not supported yet.
            # field_names = [name for name in self._scanner.field_names()]
            # if columns is not None:
            #     scan_kwargs["fields"] = [col for col in columns if col in field_names]

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

    def __init__(
        self,
        source: str | pathlib.Path | Callable[[], IO[Any]],
        bed_schema: str = "bed3+",
        compressed: bool = False,
        *,
        fields: list[str] | None = None,
        regions: str | list[str] | None = None,
        index: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        super().__init__(source, index, batch_size)

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._scanner_kwargs = dict(bed_schema=bed_schema, compressed=compressed)
        self._schema_kwargs = dict(fields=fields)

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._src,
            regions=regions,
            index=self._index_src,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
            **self._schema_kwargs,
        )


def from_bed(
    source: str | pathlib.Path | Callable[[], IO[Any]],
    bed_schema: str = "bed3+",
    compressed: bool = False,
    *,
    fields: list[str] | None = None,
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[Any]] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> BedFile:
    """
    Create a BED file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the BED file, or a callable that opens the file
        as a file-like object.
    bed_schema : str, optional
        Schema for the BED file format, by default "bed3+".
    compressed : bool, optional
        Whether the source is compressed, by default False.
    fields : list[str], optional
        Names of the fields to project.
    regions : list[str], optional
        Genomic regions to query.
    index : str, pathlib.Path, or Callable, optional
        Index file for the BED file, by default None.
    batch_size : int, optional
        Size of the batch to read.

    Returns
    -------
    BedFile

    See also
    --------
    from_gtf : Create a GTF file data source.
    from_gff : Create a GFF file data source.
    from_bigbed : Create a BigBed file data source.
    from_bigwig : Create a BigWig file data source.
    """
    return BedFile(
        source=source,
        bed_schema=bed_schema,
        compressed=compressed,
        fields=fields,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )
