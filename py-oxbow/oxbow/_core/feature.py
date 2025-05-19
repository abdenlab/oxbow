"""
This module defines classes for working with genomic feature files, including BED, BigBed, BigWig, GFF, and GTF formats.

Classes
-------
FeatureFile
    Base class for genomic feature files.
BbiFile
    Base class for BigBed and BigWig files.
BigBedFile
    Class for handling BigBed files.
BigWigFile
    Class for handling BigWig files.
BedFile
    Class for handling BED files.
GxfFile
    Base class for GFF and GTF files.
GffFile
    Class for handling GFF files.
GtfFile
    Class for handling GTF files.

Functions
---------
from_bed(uri, opener, fields, index, regions, bed_schema, compressed)
    Create a BedFile instance.
from_bigbed(uri, opener, fields, regions, resolution, bed_schema)
    Create a BigBedFile instance.
from_bigwig(uri, opener, fields, regions, resolution)
    Create a BigWigFile instance.
from_gff(uri, opener, fields, index, compressed, regions, attribute_defs, attribute_scan_rows)
    Create a GffFile instance.
from_gtf(uri, opener, fields, index, compressed, regions, attribute_defs, attribute_scan_rows)
    Create a GtfFile instance.
"""

from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any, Callable, Generator

import pyarrow as pa

from oxbow._filetypes import FileType
from oxbow._core.base import DataSource
from oxbow.oxbow import (
    PyBBIZoomScanner,
    PyBigBedScanner,
    PyBigWigScanner,
    PyGffScanner,
    PyGtfScanner,
)


class FeatureFile(DataSource):
    """
    Base class for genomic feature files.

    This class provides common functionality for handling various genomic file formats.

    Functions
    ---------
    from_bed(uri, opener, fields, index, regions, bed_schema, compressed)
        Create a BedFile instance.
    from_bigbed(uri, opener, fields, regions, resolution, bed_schema)
        Create a BigBedFile instance.
    from_bigwig(uri, opener, fields, regions, resolution)
        Create a BigWigFile instance.
    from_gff(uri, opener, fields, index, compressed, regions, attribute_defs, attribute_scan_rows)
        Create a GffFile instance.
    from_gtf(uri, opener, fields, index, compressed, regions, attribute_defs, attribute_scan_rows)
        Create a GtfFile instance.
    """

    if TYPE_CHECKING:
        from_bed: BedFile
        from_bigbed: BigBedFile
        from_bigwig: BigWigFile
        from_gff: GffFile
        from_gtf: GtfFile


class BbiFile(FeatureFile):
    """
    Base class for BigBed and BigWig files.

    This class provides functionality for handling BigBed and BigWig files, including zoom levels and region-based queries.

    Attributes
    ----------
    zoom_levels : list[int]
        List of available zoom levels for the file.
    """

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
                yield self._make_batch_reader(stream)
        else:
            stream = self._scanner.scan
            yield self._make_batch_reader(stream)

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
    """
    Class for handling BigBed files.

    Parameters
    ----------
    bed_schema : str, optional
        Schema for the BED file format, by default "bed3+".
    """

    if TYPE_CHECKING:
        _scanner: FileType.BigBed.value

    def __init__(self, *args, bed_schema="bed3+", **kwargs):
        self.__scanner_kwargs = dict(schema=bed_schema)
        super().__init__(*args, **kwargs)


class BigWigFile(BbiFile, file_type=FileType.BigWig):
    """
    Class for handling BigWig files.
    """

    if TYPE_CHECKING:
        _scanner: FileType.BigWig.value


class BedFile(FeatureFile, file_type=FileType.BED):
    """
    Class for handling BED files.

    Parameters
    ----------
    uri : str
        The URI of the BED file.
    opener : Callable, optional
        A callable to open the file.
    fields : list[str], optional
        Names of the fields to project.
    index : str, optional
        Path to the index file.
    regions : list[str], optional
        Genomic regions to query.
    bed_schema : str, optional
        Schema for the BED file format, by default "bed3".
    compressed : bool, optional
        Whether the source is compressed, by default False.
    """

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
                yield self._make_batch_reader(stream)
        else:
            stream = self._scanner.scan
            yield self._make_batch_reader(stream)

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


class GxfFile(FeatureFile):
    """
    Base class for GFF and GTF files.

    Functions
    ----------
    from_gff(uri, opener, fields, index, compressed, regions, attribute_defs, attribute_scan_rows)
        Create a GffFile instance.
    from_gtf(uri, opener, fields, index, compressed, regions, attribute_defs, attribute_scan_rows)
        Create a GtfFile instance.
    """

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
                yield self._make_batch_reader(stream)
        else:
            stream = partial(
                self._scanner.scan,
                **self._scan_kwargs,
            )
            yield self._make_batch_reader(stream)

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
        self.__scanner_kwargs = dict(compressed=compressed)
        super().__init__(uri, opener, fields)
        if attribute_defs is None:
            attribute_defs = self._scanner.attribute_defs(attribute_scan_rows)
        self.__scan_kwargs = dict(attribute_defs=attribute_defs, fields=fields)
        self.__schema_kwargs = dict(fields=fields, attribute_defs=attribute_defs)


class GffFile(GxfFile, file_type=FileType.GFF):
    """
    Class for handling GFF files.
    """

    if TYPE_CHECKING:
        _scanner: FileType.GFF.value


class GtfFile(GxfFile, file_type=FileType.GTF):
    """
    Class for handling GTF files.
    """

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
    uri : str
        The URI of the BED file.
    opener : Callable, optional
        A callable to open the file.
    fields : list[str], optional
        Names of the fields to project.
    index : str, optional
        Path to the index file.
    regions : list[str], optional
        Genomic regions to query.
    bed_schema : str, optional
        Schema for the BED file format, by default "bed3".
    compressed : bool, optional
        Whether the source is compressed, by default False.

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
    uri : str
        The URI of the BigBed file.
    opener : Callable, optional
        A callable to open the file.
    fields : list[str], optional
        Names of the fields to project.
    regions : list[str], optional
        Genomic regions to query.
    resolution : int, optional
        Resolution level for zoomed data.
    bed_schema : str, optional
        Schema for the BED file format, by default "bed3+".

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
    uri : str
        The URI of the BigWig file.
    opener : Callable, optional
        A callable to open the file.
    fields : list[str], optional
        Names of the fields to project.
    regions : list[str], optional
        Genomic regions to query.
    resolution : int, optional
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
    uri : str
        The URI of the GFF file.
    opener : Callable, optional
        A callable to open the file.
    fields : list[str], optional
        Names of the fields to project.
    index : str, optional
        Path to the index file.
    compressed : bool, optional
        Whether the source is compressed, by default False.
    regions : list[str], optional
        Genomic regions to query.
    attribute_defs : dict, optional
        Attribute definitions for the file.
    attribute_scan_rows : int, optional
        Number of rows to scan for attribute definitions, by default 1024.

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
    uri : str
        The URI of the GTF file.
    opener : Callable, optional
        A callable to open the file.
    fields : list[str], optional
        Names of the fields to project.
    index : str, optional
        Path to the index file.
    compressed : bool, optional
        Whether the source is compressed, by default False.
    regions : list[str], optional
        Genomic regions to query.
    attribute_defs : dict, optional
        Attribute definitions for the file.
    attribute_scan_rows : int, optional
        Number of rows to scan for attribute definitions, by default 1024.

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
