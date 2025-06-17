"""
DataSource classes for htslib variant call formats.
"""

from __future__ import annotations

import pathlib
from typing import IO, Callable, Generator, Literal

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

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
        source: str | Callable[[], IO[bytes] | str],
        compressed: bool = False,
        *,
        fields=None,
        info_fields: list[str] | None = None,
        samples: list[str] | None = None,
        genotype_fields: list[str] | None = None,
        genotype_by: Literal["sample", "field"] = "sample",
        regions: str | list[str] | None = None,
        index: str | Callable[[], IO[bytes] | str] | None = None,
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
        """List of INFO field definitions declared in the header."""
        return self.scanner().info_field_defs()

    @property
    def genotype_field_defs(self) -> list[tuple[str, str, str]]:
        """List of FORMAT field definitions declared in the header."""
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
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    compression: Literal["infer", "bgzf", "gzip", None] = "infer",
    *,
    fields: list[str] | None = None,
    info_fields: list[str] | None = None,
    samples: list[str] | None = None,
    genotype_fields: list[str] | None = None,
    genotype_by: Literal["sample", "field"] = "sample",
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> VcfFile:
    """
    Create a VCF file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the VCF file, or a callable that opens the file
        as a file-like object.
    compression : Literal["infer", "bgzf", "gzip", None], default: "infer"
        Compression of the source bytestream. If "infer" and ``source`` is a
        URI or path, the file's compression is guessed based on the extension,
        where ".gz" or ".bgz" is interpreted as BGZF. Pass "gzip" to decode
        regular GZIP. If None, the source bytestream is assumed to be
        uncompressed. For more customized decoding, provide a callable
        ``source`` instead.
    fields : list[str], optional
        Specific fixed fields to project. By default, all fixed fields are
        included.
    info_fields : list[str], optional [default: None]
        INFO fields to project. These will be nested under an "info" column.
        If None, all INFO fields declared in the header are included. To omit
        all INFO fields, set ``info_fields=[]``.
    samples : list[str], optional [default: None]
        A subset of samples to include in the genotype output. If None, all
        samples declared in the header are included. To omit all sample
        genotype data, set ``samples=[]``.
    genotype_fields : list[str], optional [default: None]
        Genotype (aka "FORMAT") fields to project for each sample. If None, all
        FORMAT fields declared in the header are included.
    genotype_by : Literal["sample", "field"], optional [default: "sample"]
        Determines how genotype-specific data is organized. If "sample", each
        sample is provided as a separate column with nested FORMAT fields. If
        "field", each FORMAT field is provided as a separate column with nested
        sample name fields.
    regions : str | list[str], optional
        One or more genomic regions to query. Only applicable if an associated
        index file is available.
    index : str, pathlib.Path, or Callable, optional
        An optional index file associated with the VCF file. If ``source`` is a
        URI or path, is BGZF-compressed, and the index file shares the same
        name with a ".tbi" or ".csi" extension, the index file is automatically
        detected.
    batch_size : int, optional [default: 131072]
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
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    compression: Literal["bgzf", None] = "bgzf",
    *,
    fields: list[str] | None = None,
    info_fields: list[str] | None = None,
    samples: list[str] | None = None,
    genotype_fields: list[str] | None = None,
    genotype_by: Literal["sample", "field"] = "sample",
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> BcfFile:
    """
    Create a BCF file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the BCF file, or a callable that opens the file
        as a file-like object.
    compression : Literal["bgzf", None], default: "bgzf"
        Compression of the source bytestream. By default, BCF sources are
        assumed to be BGZF-compressed. If None, the source is assumed to be
        uncompressed. For more custom decoding, provide a callable ``source``
        instead.
    fields : list[str], optional
        Specific fixed fields to project. By default, all fixed fields are
        included.
    info_fields : list[str], optional [default: None]
        INFO fields to project. These will be nested under an "info" column.
        If None, all INFO fields declared in the header are included. To omit
        all INFO fields, set ``info_fields=[]``.
    samples : list[str], optional [default: None]
        A subset of samples to include in the genotype output. If None, all
        samples declared in the header are included. To omit all sample
        genotype data, set ``samples=[]``.
    genotype_fields : list[str], optional [default: None]
        Genotype (aka "FORMAT") fields to project for each sample. If None, all
        FORMAT fields declared in the header are included.
    genotype_by : Literal["sample", "field"], optional [default: "sample"]
        Determines how genotype-specific data is organized. If "sample", each
        sample is provided as a separate column with nested FORMAT fields. If
        "field", each FORMAT field is provided as a separate column with nested
        sample name fields.
    regions : str | list[str], optional
        One or more genomic regions to query. Only applicable if an associated
        index file is available.
    index : str, optional
        An optional index file associated with the BCF file. If ``source`` is a
        URI or path, is BGZF-compressed, and the index file shares the same
        name with a ".csi" extension, the index file is automatically
        detected.
    batch_size : int, optional [default: 131072]
        The number of records to read in each batch.

    Returns
    -------
    BcfFile
        A data source object representing the BCF file.

    Notes
    -----
    The Binary Call Format (BCF) is a binary representation of the Variant Call
    Format (VCF), designed for efficient storage and processing of genomic
    variant data. It is commonly used in large-scale sequencing projects.

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
