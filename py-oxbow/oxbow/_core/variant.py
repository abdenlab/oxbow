"""
DataSource classes for htslib variant call formats.
"""

from __future__ import annotations

import pathlib
from typing import IO, Any, Callable, Generator, Literal, Self

import pyarrow as pa

from oxbow._core.base import DEFAULT_BATCH_SIZE, DataSource, prepare_source_and_index
from oxbow.oxbow import PyBcfScanner, PyVcfScanner


class VariantFile(DataSource):
    def _batchreader_builder(
        self,
        scan_fn,
        field_names: list[str],
        region: str | None = None,
    ) -> Callable[[list[str] | None, int], pa.RecordBatchReader]:
        def builder(columns, batch_size):
            scan_kwargs = self._schema_kwargs.copy()

            if columns is not None:
                n = len(columns)
                scan_kwargs["fields"] = [col for col in columns if col in field_names]
                n -= len(scan_kwargs["fields"])

                if "info" not in columns:
                    scan_kwargs["info_fields"] = []
                elif scan_kwargs.get("info_fields") == []:
                    raise ValueError(
                        "Cannot select `info` column if no info fields are provided."
                    )
                else:
                    n -= 1
                if n == 0:
                    scan_kwargs["samples"] = []

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
        source: str | Callable[[], IO[Any] | str],
        compressed: bool = False,
        *,
        fields=None,
        info_fields: list[str] | None = None,
        samples: list[str] | None = None,
        genotype_fields: list[str] | None = None,
        genotype_by: Literal["sample", "field"] = "sample",
        regions: str | list[str] | None = None,
        index: str | Callable[[], IO[Any] | str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        super().__init__(source, index, batch_size)

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._scanner_kwargs = dict(compressed=compressed)
        self._schema_kwargs = dict(
            fields=fields,
            info_fields=info_fields,
            samples=samples,
            genotype_fields=genotype_fields,
            genotype_by=genotype_by,
        )

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._src,
            regions=regions,
            index=self._index_src,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
            **self._schema_kwargs,
        )

    @property
    def chrom_names(self) -> list[str]:
        """List of reference sequence names declared in the header."""
        return self.scanner().chrom_names()

    @property
    def chrom_sizes(self) -> list[tuple[str, int]]:
        """List of reference sequence names and their lengths in bp."""
        return self.scanner().chrom_sizes()

    @property
    def info_field_defs(self) -> list[tuple[str, str, str]]:
        """List of INFO field definitions in the VCF header."""
        return self.scanner().info_field_defs()

    @property
    def genotype_field_defs(self) -> list[tuple[str, str, str]]:
        """List of FORMAT field definitions in the VCF header."""
        return self.scanner().genotype_field_defs()

    @property
    def samples(self) -> list[str]:
        """List of sample IDs declared in the header."""
        return self.scanner().sample_names()


class VcfFile(VariantFile):
    _scanner_type = PyVcfScanner


class BcfFile(VariantFile):
    _scanner_type = PyBcfScanner


def from_vcf(
    source: str | pathlib.Path | Callable[[], IO[Any] | str],
    compression: Literal["infer", "bgzf", "gzip", None] = "infer",
    *,
    fields: list[str] | None = None,
    info_fields: list[str] | None = None,
    samples: list[str] | None = None,
    genotype_fields: list[str] | None = None,
    genotype_by: Literal["sample", "field"] = "sample",
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[Any] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> VcfFile:
    """
    Create a VCF (Variant Call Format) file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the VCF file, or a callable that opens the file
        as a file-like object.
    compression : Literal["infer", "bgzf", "gzip", None], default: "infer"
        If "infer" and `source` is a URI or path, the file's compression is
        guessed based on the file extension, where ".gz" or ".bgz" is
        interpreted as BGZF. To decode vanilla GZIP, use "gzip". If None, the
        source bytestream is assumed to be uncompressed. For more custom
        decoding, provide a callable `source` instead.
    fields : list[str], optional
        Specific fields to include from the VCF file.
    info_fields : list[str], optional
        INFO fields to extract from the VCF file.
    samples : list[str], optional
        A subset of samples to include.
    genotype_fields : list[str], optional
        FORMAT fields to extract genotype-specific information.
    genotype_by : Literal["sample", "field"], optional
        Determines how genotype data is organized, by default "sample".
    regions : list[str], optional
        Genomic regions to query, specified as strings (e.g., "chr1:1000-2000").
    index : str, pathlib.Path, or Callable, optional
        The index file associated with the VCF file.
    batch_size : int, optional
        The number of records to read in each batch.

    Returns
    -------
    VcfFile
        A data source object representing the VCF file.

    Notes
    -----
    The Variant Call Format (VCF) is a text-based format used to store
    information about genomic variants. It is widely used in bioinformatics
    for storing and sharing variant data from sequencing projects.

    See Also
    --------
    from_bcf : Create a BCF file data source.
    """
    source, index, bgzf_compressed = prepare_source_and_index(
        source, index, compression
    )
    return VcfFile(
        source=source,
        compressed=bgzf_compressed,
        fields=fields,
        info_fields=info_fields,
        samples=samples,
        genotype_fields=genotype_fields,
        genotype_by=genotype_by,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )


def from_bcf(
    source: str | pathlib.Path | Callable[[], IO[Any] | str],
    compression: Literal["bgzf", None] = "bgzf",
    *,
    fields: list[str] | None = None,
    info_fields: list[str] | None = None,
    samples: list[str] | None = None,
    genotype_fields: list[str] | None = None,
    genotype_by: Literal["sample", "field"] = "sample",
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[Any] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> BcfFile:
    """
    Create a BCF (Binary Call Format) file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the VCF file, or a callable that opens the file
        as a file-like object.
    compression : Literal["bgzf", None], default: "bgzf"
        Compression of the source bytestream. By default, BCF sources are
        assumed to be BGZF-compressed. If None, the source is assumed to be
        uncompressed. For more custom decoding, provide a callable `source`
        instead.
    fields : list[str], optional
        Specific fields to include from the BCF file.
    info_fields : list[str], optional
        INFO fields to extract from the BCF file.
    samples : list[str], optional
        A subset of samples to include.
    genotype_fields : list[str], optional
        FORMAT fields to extract genotype-specific information.
    genotype_by : Literal["sample", "field"], optional
        Determines how genotype data is organized, by default "sample".
    regions : list[str], optional
        Genomic regions to query, specified as strings (e.g., "chr1:1000-2000").
    index : str, optional
        The index file associated with the BCF file.
    batch_size : int, optional
        The number of records to read in each batch.

    Returns
    -------
    BcfFile
        A data source object representing the BCF file.

    Notes
    -----
    The BCF format is a binary representation of the Variant Call Format (VCF),
    designed for efficient storage and processing of genomic variant data.
    It is commonly used in large-scale sequencing projects.

    See Also
    --------
    from_vcf : Create a VCF file data source.
    """
    source, index, bgzf_compressed = prepare_source_and_index(
        source, index, compression
    )
    return BcfFile(
        source=source,
        compressed=bgzf_compressed,
        fields=fields,
        info_fields=info_fields,
        samples=samples,
        genotype_fields=genotype_fields,
        genotype_by=genotype_by,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )
