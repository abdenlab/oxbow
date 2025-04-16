from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any, Callable, Generator

import pyarrow as pa

from oxbow._filetypes import FileType
from oxbow._core.base import DataFile
from oxbow.oxbow import (
    PyBigBedScanner,
    PyBigWigScanner,
    PyBBIZoomScanner,
    PyGffScanner,
    PyGtfScanner,
)


class FunctionFile(DataFile):
    if TYPE_CHECKING:
        from_bed: BedFile
        from_bigbed: BigBedFile
        from_bigwig: BigWigFile
        from_gff: GffFile
        from_gtf: GtfFile


class BbiFile(FunctionFile):
    if TYPE_CHECKING:
        from_bigbed: BigBedFile
        from_bigwig: BigWigFile

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
                stream = partial(self._scanner.scan_query, region=region)
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

    @property
    def _scanner(
        self,
    ) -> PyBigBedScanner | PyBigWigScanner | PyBBIZoomScanner:
        scanner = super()._scanner
        assert isinstance(scanner, (PyBigBedScanner, PyBigWigScanner))
        if self._resolution:
            return scanner.get_zoom(self._resolution)
        else:
            return scanner

    def __init__(
        self,
        uri,
        opener=None,
        fields=None,
        regions=None,
        resolution=None,
    ):
        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions
        self._resolution = resolution
        self.__schema_kwargs = dict(fields=fields)
        self.__scan_kwargs = dict(fields=fields)
        self.__scanner_kwargs = {}
        super().__init__(uri, opener, fields)

    @property
    def zoom_levels(self):
        scanner = super()._scanner
        assert isinstance(scanner, (PyBigBedScanner, PyBigWigScanner))
        return scanner.zoom_levels()

    def zoom(self, resolution: int) -> BbiFile:
        """
        Create a new BbiFile instance with the specified resolution.

        Parameters
        ----------
        resolution
            The resolution level for zoomed data.

        Returns
        -------
        BbiFile
            A new instance of the BbiFile class with the specified resolution.
        """
        return type(self)(
            self._uri,
            self._opener,
            **dict(
                **self._scan_kwargs,
                **self._scanner_kwargs,
                **self._schema_kwargs,
                resolution=resolution,
            ),
        )


class BigBedFile(BbiFile, file_type=FileType.BigBed):
    if TYPE_CHECKING:
        _scanner: FileType.BigBed.value

    def __init__(self, *args, bed_schema="bed3+", **kwargs):
        self.__scanner_kwargs = dict(schema=bed_schema)
        super().__init__(*args, **kwargs)


class BigWigFile(BbiFile, file_type=FileType.BigWig):
    if TYPE_CHECKING:
        _scanner: FileType.BigWig.value


class BedFile(FunctionFile, file_type=FileType.BED):
    if TYPE_CHECKING:
        _scanner: FileType.BED.value

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
                stream = partial(
                    self._scanner.scan_query, region=region, index=self._index
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
        regions=None,
        bed_schema="bed3",
        compressed=False,
    ):
        self._index = index
        self._regions = regions
        self.__schema_kwargs = dict(fields=fields)
        self.__scan_kwargs = dict(fields=fields)
        self.__scanner_kwargs = dict(bed_schema=bed_schema, compressed=compressed)
        super().__init__(uri, opener, fields)


class GxfFile(FunctionFile):
    if TYPE_CHECKING:
        from_gff: GffFile
        from_gtf: GtfFile
        _scanner: PyGffScanner | PyGtfScanner

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
        uri,
        opener=None,
        fields=None,
        index=None,
        compressed=False,
        regions=None,
        attribute_defs=None,
        attribute_scan_rows=1024,
    ):
        self._index = index
        self._regions = regions
        self.__schema_kwargs = dict(fields=fields)
        self.__scanner_kwargs = dict(compressed=compressed)
        super().__init__(uri, opener, fields)
        if attribute_defs is None:
            attribute_defs = self._scanner.attribute_defs(attribute_scan_rows)
        self.__scan_kwargs = dict(attribute_defs=attribute_defs, fields=fields)


class GffFile(GxfFile, file_type=FileType.GFF):
    if TYPE_CHECKING:
        _scanner: FileType.GFF.value


class GtfFile(GxfFile, file_type=FileType.GTF):
    if TYPE_CHECKING:
        _scanner: FileType.GTF.value


def from_bed(
    uri: str,
    opener: Callable | None = None,
    fields: list[str] | None = None,
    index: str | None = None,
    regions: list[str] | None = None,
    bed_schema: str = "bed3",
    compressed: bool = False,
) -> BedFile:
    """
    Create a BedFile instance.

    Parameters
    ----------
    uri
        The URI of the BED file.
    opener
        A callable to open the file.
    fields
        Names of the fields to project.
    index
        Path to the index file.
    regions
        Genomic regions to query.
    bed_schema
        Schema for the BED file format.
    compressed
        Whether the source is compressed.

    Returns
    -------
    BedFile
        An instance of the BedFile class.
    """
    return BedFile(
        uri,
        opener=opener,
        fields=fields,
        index=index,
        regions=regions,
        bed_schema=bed_schema,
        compressed=compressed,
    )


def from_bigbed(
    uri: str,
    opener: Callable | None = None,
    fields: list[str] | None = None,
    regions: list[str] | None = None,
    resolution: int | None = None,
    bed_schema: str = "bed3+",
) -> BigBedFile:
    """
    Create a BigBedFile instance.

    Parameters
    ----------
    uri
        The URI of the BigBed file.
    opener
        A callable to open the file.
    fields
        Names of the fields to project.
    regions
        Genomic regions to query.
    resolution
        Resolution level for zoomed data.
    bed_schema
        Schema for the BED file format.

    Returns
    -------
    BigBedFile
        An instance of the BigBedFile class.
    """
    return BigBedFile(
        uri,
        opener=opener,
        fields=fields,
        regions=regions,
        resolution=resolution,
        bed_schema=bed_schema,
    )


def from_bigwig(
    uri: str,
    opener: Callable | None = None,
    fields: list[str] | None = None,
    regions: list[str] | None = None,
    resolution: int | None = None,
) -> BigWigFile:
    """
    Create a BigWigFile instance.

    Parameters
    ----------
    uri
        The URI of the BigWig file.
    opener
        A callable to open the file.
    fields
        Names of the fields to project.
    regions
        Genomic regions to query.
    resolution
        Resolution level for zoomed data.

    Returns
    -------
    BigWigFile
        An instance of the BigWigFile class.
    """
    return BigWigFile(
        uri,
        opener=opener,
        fields=fields,
        regions=regions,
        resolution=resolution,
    )


def from_gff(
    uri: str,
    opener: Callable | None = None,
    fields: list[str] | None = None,
    index: str | None = None,
    compressed: bool = False,
    regions: list[str] | None = None,
    attribute_defs: dict | None = None,
    attribute_scan_rows: int = 1024,
) -> GffFile:
    """
    Create a GffFile instance.

    Parameters
    ----------
    uri
        The URI of the GFF file.
    opener
        A callable to open the file.
    fields
        Names of the fields to project.
    index
        Path to the index file.
    compressed
        Whether the source is compressed.
    regions
        Genomic regions to query.
    attribute_defs
        Attribute definitions for the file.
    attribute_scan_rows
        Number of rows to scan for attribute definitions.

    Returns
    -------
    GffFile
        An instance of the GffFile class.
    """
    return GffFile(
        uri,
        opener=opener,
        fields=fields,
        index=index,
        compressed=compressed,
        regions=regions,
        attribute_defs=attribute_defs,
        attribute_scan_rows=attribute_scan_rows,
    )


def from_gtf(
    uri: str,
    opener: Callable | None = None,
    fields: list[str] | None = None,
    index: str | None = None,
    compressed: bool = False,
    regions: list[str] | None = None,
    attribute_defs: dict | None = None,
    attribute_scan_rows: int = 1024,
) -> GtfFile:
    """
    Create a GtfFile instance.

    Parameters
    ----------
    uri
        The URI of the GTF file.
    opener
        A callable to open the file.
    fields
        Names of the fields to project.
    index
        Path to the index file.
    compressed
        Whether the source is compressed.
    regions
        Genomic regions to query.
    attribute_defs
        Attribute definitions for the file.
    attribute_scan_rows
        Number of rows to scan for attribute definitions.

    Returns
    -------
    GtfFile
        An instance of the GtfFile class.
    """
    return GtfFile(
        uri,
        opener=opener,
        fields=fields,
        index=index,
        compressed=compressed,
        regions=regions,
        attribute_defs=attribute_defs,
        attribute_scan_rows=attribute_scan_rows,
    )
