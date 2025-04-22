"""
This module defines classes and functions for working with variant files, including VCF and BCF formats.

Classes
-------
VariantFile
    Base class for variant files.
BcfFile
    Class for handling BCF files.
VcfFile
    Class for handling VCF files.

Functions
---------
from_bcf(uri, opener, fields, index, regions, info_fields, genotype_fields, samples, genotype_by)
    Create a BcfFile instance from a BCF file.
from_vcf(uri, opener, fields, index, regions, info_fields, genotype_fields, samples, genotype_by)
    Create a VcfFile instance from a VCF file.
"""

from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any, Callable, Generator, Literal

import pyarrow as pa

from oxbow._core.base import DataFile
from oxbow._filetypes import FileType
from oxbow.oxbow import PyBcfScanner, PyVcfScanner


class VariantFile(DataFile):
    """
    Base class for variant files.

    This class provides common functionality for handling VCF and BCF files.

    Functions
    ---------
    from_bcf(uri, opener, fields, index, regions, info_fields, genotype_fields, samples, genotype_by)
        Create a BcfFile instance from a BCF file.
    from_vcf(uri, opener, fields, index, regions, info_fields, genotype_fields, samples, genotype_by)
        Create a VcfFile instance from a VCF file.
    """

    if TYPE_CHECKING:
        from_bcf: BcfFile
        from_vcf: VcfFile
        _scanner: PyBcfScanner | PyVcfScanner

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
            stream = partial(self._scanner.scan, **self._scan_kwargs)
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
        info_fields: list[str] | None = None,
        genotype_fields: list[str] | None = None,
        samples: list[str] | None = None,
        genotype_by: Literal["sample", "field"] = "sample",
    ):
        self._index = index
        self._regions = regions
        super().__init__(uri, opener, fields)
        self.__scanner_kwargs = dict(compressed=compressed)
        if info_fields is None:
            info_fields = self._scanner.info_field_names()
        else:
            info_fields = [
                d for d in self._scanner.info_field_names() if d in info_fields
            ]
        if genotype_fields is None:
            genotype_fields = self._scanner.genotype_field_names()
        else:
            genotype_fields = [
                d for d in self._scanner.genotype_field_names() if d in genotype_fields
            ]
        if samples is None:
            samples = self._scanner.sample_names()
        else:
            samples = [s for s in self._scanner.sample_names() if s in samples]
        self.__schema_kwargs = dict(
            fields=fields,
            genotype_by=genotype_by,
            genotype_fields=genotype_fields,
            info_fields=info_fields,
            samples=samples,
        )
        self.__scan_kwargs = dict(
            fields=fields,
            genotype_by=genotype_by,
            genotype_fields=genotype_fields,
            info_fields=info_fields,
            samples=samples,
        )


class BcfFile(VariantFile, file_type=FileType.BCF):
    """
    Class for handling BCF files.

    Parameters
    ----------
    uri : str
        The URI or file path to the BCF file.
    opener : Callable, optional
        A custom file opener, such as one for handling remote files.
    fields : list[str], optional
        Specific fields to include from the BCF file.
    index : str, optional
        Path to the index file associated with the BCF file.
    regions : list[str], optional
        Genomic regions to query, specified as strings (e.g., "chr1:1000-2000").
    info_fields : list[str], optional
        INFO fields to extract from the BCF file.
    genotype_fields : list[str], optional
        FORMAT fields to extract genotype-specific information.
    samples : list[str], optional
        A subset of samples to include.
    genotype_by : Literal["sample", "field"], optional
        Determines how genotype data is organized, by default "sample".
    """

    if TYPE_CHECKING:
        _scanner: FileType.BCF.value

    def __init__(
        self,
        uri,
        opener=None,
        fields=None,
        index=None,
        regions=None,
        info_fields: list[str] | None = None,
        genotype_fields: list[str] | None = None,
        samples: list[str] | None = None,
        genotype_by: Literal["sample"] | Literal["field"] = "sample",
    ):
        super().__init__(
            uri=uri,
            opener=opener,
            fields=fields,
            index=index,
            compressed=True,
            regions=regions,
            info_fields=info_fields,
            genotype_fields=genotype_fields,
            samples=samples,
            genotype_by=genotype_by,
        )

    @property
    def _scanner_kwargs(self) -> dict[str, Any]:
        return {}


class VcfFile(VariantFile, file_type=FileType.VCF):
    """
    Class for handling VCF files.

    Parameters
    ----------
    uri : str
        The URI or file path to the VCF file.
    opener : Callable, optional
        A custom file opener, such as one for handling remote files.
    fields : list[str], optional
        Specific fields to include from the VCF file.
    index : str, optional
        Path to the index file associated with the VCF file.
    regions : list[str], optional
        Genomic regions to query, specified as strings (e.g., "chr1:1000-2000").
    info_fields : list[str], optional
        INFO fields to extract from the VCF file.
    genotype_fields : list[str], optional
        FORMAT fields to extract genotype-specific information.
    samples : list[str], optional
        A subset of samples to include.
    genotype_by : Literal["sample", "field"], optional
        Determines how genotype data is organized, by default "sample".
    """

    if TYPE_CHECKING:
        _scanner: FileType.VCF.value


def from_bcf(
    uri: str,
    opener: Callable | None = None,
    fields: list[str] | None = None,
    index: str | None = None,
    regions: list[str] | None = None,
    info_fields: list[str] | None = None,
    genotype_fields: list[str] | None = None,
    samples: list[str] | None = None,
    genotype_by: Literal["sample", "field"] = "sample",
) -> BcfFile:
    """
    Create a `BcfFile` object from a BCF (Binary Call Format) file.

    Parameters
    ----------
    uri : str
        The URI or file path to the BCF file.
    opener : Callable, optional
        A custom file opener, such as one for handling remote files.
    fields : list[str], optional
        Specific fields to include from the BCF file.
    index : str, optional
        Path to the index file associated with the BCF file.
    regions : list[str], optional
        Genomic regions to query, specified as strings (e.g., "chr1:1000-2000").
    info_fields : list[str], optional
        INFO fields to extract from the BCF file.
    genotype_fields : list[str], optional
        FORMAT fields to extract genotype-specific information.
    samples : list[str], optional
        A subset of samples to include.
    genotype_by : Literal["sample", "field"], optional
        Determines how genotype data is organized, by default "sample".

    Returns
    -------
    BcfFile
        An object representing the BCF file, ready for querying and analysis.

    Notes
    -----
    The BCF format is a binary representation of the Variant Call Format (VCF), designed
    for efficient storage and processing of genomic variant data. It is commonly used in
    large-scale sequencing projects.

    See Also
    --------
    from_vcf : Create a `VcfFile` object from a VCF file.
    """
    return BcfFile(
        uri=uri,
        opener=opener,
        fields=fields,
        index=index,
        regions=regions,
        info_fields=info_fields,
        genotype_fields=genotype_fields,
        samples=samples,
        genotype_by=genotype_by,
    )


def from_vcf(
    uri: str,
    opener: Callable | None = None,
    fields: list[str] | None = None,
    index: str | None = None,
    regions: list[str] | None = None,
    info_fields: list[str] | None = None,
    genotype_fields: list[str] | None = None,
    samples: list[str] | None = None,
    genotype_by: Literal["sample", "field"] = "sample",
) -> VcfFile:
    """
    Create a `VcfFile` object from a VCF (Variant Call Format) file.

    Parameters
    ----------
    uri : str
        The URI or file path to the VCF file.
    opener : Callable, optional
        A custom file opener, such as one for handling remote files.
    fields : list[str], optional
        Specific fields to include from the VCF file.
    index : str, optional
        Path to the index file associated with the VCF file.
    regions : list[str], optional
        Genomic regions to query, specified as strings (e.g., "chr1:1000-2000").
    info_fields : list[str], optional
        INFO fields to extract from the VCF file.
    genotype_fields : list[str], optional
        FORMAT fields to extract genotype-specific information.
    samples : list[str], optional
        A subset of samples to include.
    genotype_by : Literal["sample", "field"], optional
        Determines how genotype data is organized, by default "sample".

    Returns
    -------
    VcfFile
        An object representing the VCF file, ready for querying and analysis.

    Notes
    -----
    The Variant Call Format (VCF) is a text-based format used to store information about
    genomic variants. It is widely used in bioinformatics for storing and sharing variant
    data from sequencing projects.

    See Also
    --------
    from_bcf : Create a `BcfFile` object from a BCF file.
    """
    return VcfFile(
        uri=uri,
        opener=opener,
        fields=fields,
        index=index,
        regions=regions,
        info_fields=info_fields,
        genotype_fields=genotype_fields,
        samples=samples,
        genotype_by=genotype_by,
    )
