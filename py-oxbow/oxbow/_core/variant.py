from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any, Callable, Generator, Literal

import pyarrow as pa

from oxbow._core.base import DataFile
from oxbow._filetypes import FileType
from oxbow.oxbow import PyBcfScanner, PyVcfScanner


class VariantFile(DataFile):
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
    uri
        The URI or file path to the BCF file.
    opener
        A custom file opener, such as one for handling remote files.
    fields
        Specific fields to include from the BCF file. If None, all fields are included.
    index
        Path to the index file associated with the BCF file.
    regions
        Genomic regions to query, specified as strings (e.g., "chr1:1000-2000").
    info_fields
        INFO fields to extract from the BCF file. INFO fields provide metadata about variants,
        such as allele frequency (AF) or total depth (DP).
    genotype_fields
        FORMAT fields to extract genotype-specific information, such as genotype (GT) or
        genotype quality (GQ).
    samples
        A subset of samples to include. If None, all samples in the file are included.
    genotype_by
        Determines how genotype data is organized. If 'sample', data is grouped by sample.
        If 'field', data is grouped by field. Default is 'sample'.

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
    uri
        The URI or file path to the VCF file.
    opener
        A custom file opener, such as one for handling remote files.
    fields
        Specific fields to include from the VCF file. If None, all fields are included.
    index
        Path to the index file associated with the VCF file.
    regions
        Genomic regions to query, specified as strings (e.g., "chr1:1000-2000").
    info_fields
        INFO fields to extract from the VCF file. INFO fields provide metadata about variants,
        such as allele frequency (AF) or total depth (DP).
    genotype_fields
        FORMAT fields to extract genotype-specific information, such as genotype (GT) or
        genotype quality (GQ).
    samples
        A subset of samples to include. If None, all samples in the file are included.
    genotype_by
        Determines how genotype data is organized. If 'sample', data is grouped by sample.
        If 'field', data is grouped by field. Default is 'sample'.

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
