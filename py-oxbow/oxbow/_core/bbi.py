"""
DataSource classes for BBI (BigWig and BigBed) formats and their zoom levels.
"""

from __future__ import annotations

import pathlib
from typing import IO, Callable, Generator

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

import pyarrow as pa

from oxbow._core.base import DEFAULT_BATCH_SIZE, DataSource, prepare_source_and_index
from oxbow.oxbow import (
    PyBBIZoomScanner,
    PyBigBedScanner,
    PyBigWigScanner,
)


class BbiFile(DataSource):
    _regions: list[str] | None

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

            if region is not None:
                scan_kwargs["region"] = region

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

    @property
    def chrom_names(self) -> list[str]:
        """List of reference sequence names."""
        return self.scanner().chrom_names()

    @property
    def chrom_sizes(self) -> list[tuple[str, int]]:
        """List of reference sequence names and their lengths in bp."""
        return self.scanner().chrom_sizes()

    @property
    def zoom_levels(self):
        """List of zoom levels available."""
        return self.scanner().zoom_levels()

    def zoom(
        self,
        resolution: int,
        *,
        fields: list[str] | None = None,
        regions: str | list[str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ) -> BbiZoom:
        """
        Create a data source for a BBI file zoom level.

        Parameters
        ----------
        resolution: int
            The resolution / reduction level for zoomed data, in bp.

        Returns
        -------
        BbiZoom
            A data source representing a BBI file zoom level.
        """
        return BbiZoom(
            self, resolution, fields=fields, regions=regions, batch_size=batch_size
        )


class BigBedFile(BbiFile):
    _scanner_type = PyBigBedScanner

    def __init__(
        self,
        source: str | Callable[[], IO[bytes] | str],
        schema: str = "bed3+",
        *,
        fields: list[str] | None = None,
        regions: str | list[str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        super().__init__(source, None, batch_size)

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._schema_kwargs = dict(fields=fields)
        self._scanner_kwargs = dict(schema=schema)

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._src,
            regions=regions,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
            **self._schema_kwargs,
        )


class BigWigFile(BbiFile):
    _scanner_type = PyBigWigScanner

    def __init__(
        self,
        source: str | Callable[[], IO[bytes] | str],
        *,
        fields: list[str] | None = None,
        regions: str | list[str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        super().__init__(source, None, batch_size)

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._schema_kwargs = dict(fields=fields)
        self._scanner_kwargs = {}

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._src,
            regions=regions,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
            **self._schema_kwargs,
        )


class BbiZoom(DataSource):
    _scanner_type = PyBBIZoomScanner

    def scanner(self) -> PyBBIZoomScanner:
        return self._base.scanner().get_zoom(self._resolution)

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

            if region is not None:
                scan_kwargs["region"] = region

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
        base: BbiFile,
        resolution: int,
        *,
        fields: list[str] | None = None,
        regions: str | list[str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        if isinstance(regions, str):
            regions = [regions]

        self._batch_size = batch_size
        self._base = base
        self._resolution = resolution
        self._regions = regions
        self._schema_kwargs = dict(fields=fields)

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._base,
            self._resolution,
            regions=regions,
            batch_size=self._batch_size,
            **self._schema_kwargs,
        )


def from_bigbed(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    schema: str = "bed3+",
    *,
    fields: list[str] | None = None,
    regions: str | list[str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> BigBedFile:
    """
    Create a BigBed file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the BigBed file, or a callable that opens the file
        as a file-like object.
    bed_schema : str, optional [default: "bed3+"]
        Schema for intepreting the BED fields. The default is "bed3+", which
        includes the first three standard fields (chrom, start, end) and
        any additional data is lumped into a single "rest" column. If the
        BigBed file contains an AutoSql definition of its fields, pass
        "autosql" to use it.
    fields : list[str], optional
        Specific fields to project as columns. By default, all available fields
        are included.
    regions : list[str], optional
        One or more genomic regions to query. Only applicable if an associated
        index file is available.
    batch_size : int, optional [default: 131072]
        The number of records to read in each batch.

    Returns
    -------
    BigBedFile
        A data source object representing the BigBed file.

    See also
    --------
    from_bed : Create a BED file data source.
    from_bigwig : Create a BigWig file data source.
    :meth:`bbi.BigBedFile.zoom` : Create a data source for a zoom level.
    """
    source, *_ = prepare_source_and_index(source, None, None)
    return BigBedFile(
        source=source,
        schema=schema,
        fields=fields,
        regions=regions,
        batch_size=batch_size,
    )


def from_bigwig(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    *,
    fields: list[str] | None = None,
    regions: str | list[str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> BigWigFile:
    """
    Create a BigWig file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the BigWig file, or a callable that opens the file
        as a file-like object.
    fields : list[str], optional
        Specific fields to project as columns. By default, all available fields
        are included.
    regions : str | list[str], optional
        One or more genomic regions to query. Only applicable if an associated
        index file is available.
    batch_size : int, optional [default: 131072]
        The number of records to read in each batch.

    Returns
    -------
    BigWigFile
        A data source object representing the BigWig file.

    See also
    --------
    from_bed : Create a BED file data source.
    from_bigbed : Create a BigBed file data source.
    :meth:`bbi.BigWigFile.zoom` : Create a data source for a zoom level.
    """
    source, *_ = prepare_source_and_index(source, None, None)
    return BigWigFile(
        source=source,
        fields=fields,
        regions=regions,
        batch_size=batch_size,
    )
