"""
DataSource classes for htslib variant call formats.
"""

from __future__ import annotations

import pathlib
from typing import IO, Callable, Literal

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

from oxbow._core.base import DEFAULT_BATCH_SIZE, DataSource, prepare_source_and_index
from oxbow.oxbow import PyBcfScanner, PyVcfScanner


class VariantFile(DataSource):
    def __init__(
        self,
        source: str | Callable[[], IO[bytes] | str],
        compressed: bool = False,
        *,
        fields: Literal["*"] | list[str] | None = "*",
        info_fields: Literal["*"] | list[str] | None = "*",
        genotype_fields: Literal["*"] | list[str] | None = "*",
        genotype_by: Literal["sample", "field"] = "sample",
        samples: Literal["*"] | list[str] | None = None,
        samples_nested: bool = False,
        coords: Literal["01", "11"] = "11",
        regions: str | list[str] | None = None,
        index: str | Callable[[], IO[bytes] | str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        super().__init__(source, index, batch_size)

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._scanner_kwargs = dict(
            compressed=compressed,
            fields=fields,
            info_fields=info_fields,
            genotype_fields=genotype_fields,
            genotype_by=genotype_by,
            samples=samples,
            samples_nested=samples_nested,
            coords=coords,
        )

    def _scan_query(self, scanner, region, columns, batch_size):
        return scanner.scan_query(
            region=region, index=self._index, columns=columns, batch_size=batch_size
        )

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._src,
            regions=regions,
            index=self._index_src,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
        )

    def with_samples(
        self,
        samples: Literal["*"] | list[str] | None = "*",
        *,
        genotype_fields: Literal["*"] | list[str] | None = "*",
        group_by: Literal["sample", "field"] = "sample",
    ) -> Self:
        """
        Return a new data source with sample genotype data nested under a
        single ``"samples"`` struct column.

        Parameters
        ----------
        samples : "*", list[str], or None, optional [default: "*"]
            Names of samples to include in the genotype output. ``"*"``
            includes all samples declared in the header. Pass a list to select
            specific samples. ``None`` omits all sample genotype data.
        genotype_fields : "*", list[str], or None, optional [default: "*"]
            Genotype (aka FORMAT) fields to project for each sample. ``"*"``
            includes all FORMAT fields declared in the header. Pass a list to
            select specific fields. ``None`` omits all genotype fields.
        group_by : Literal["sample", "field"], optional [default: "sample"]
            Determines how genotype data is organized within the ``"samples"``
            struct. If ``"sample"``, each sample name is a sub-column with
            nested genotype fields. If ``"field"``, each genotype field is a
            sub-column with nested sample values.

        Returns
        -------
        Self
            A new data source with sample genotype data nested under a single
            ``"samples"`` struct column.
        """
        self._scanner_kwargs.update(
            genotype_fields=genotype_fields,
            genotype_by=group_by,
            samples=samples,
            samples_nested=True,
        )
        return type(self)(
            self._src,
            regions=self._regions,
            index=self._index_src,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
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
    fields: Literal["*"] | list[str] | None = "*",
    info_fields: Literal["*"] | list[str] | None = "*",
    genotype_fields: Literal["*"] | list[str] | None = "*",
    genotype_by: Literal["sample", "field"] = "sample",
    samples: Literal["*"] | list[str] | None = None,
    samples_nested: bool = False,
    coords: Literal["01", "11"] = "11",
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> VcfFile:
    """
    Create a VCF file data source.

    .. versionchanged:: 0.7.0
        The ``samples`` parameter now defaults to omitting sample genotype
        data (``None``) instead of including all samples (``"*"``). To include
        samples, pass a value to the ``samples`` parameter or use the
        :meth:`~oxbow.core.VcfFile.with_samples()` method on the returned data
        source.

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
    fields : ``"*"``, list[str], or None, optional [default: ``"*"``]
        Fixed fields to project. ``"*"`` includes all standard fields. Pass a
        list to select specific fields. ``None`` omits all fixed fields.
    info_fields : ``"*"``, list[str], or None, optional [default: ``"*"``]
        INFO fields to project, nested under an ``"info"`` column. ``"*"``
        includes all INFO fields declared in the header. Pass a list to select
        specific fields. ``None`` omits the info column entirely.
    genotype_fields : ``"*"``, list[str], or None, optional [default: ``"*"``]
        Genotype (aka "FORMAT") fields to project for each sample. ``"*"``
        includes all FORMAT fields declared in the header. Pass a list to select
        specific fields. ``None`` omits the genotype fields.
    genotype_by : Literal["sample", "field"], optional [default: "sample"]
        Determines how genotype-specific data is organized. If "sample", each
        sample is provided as a separate column with nested FORMAT fields. If
        "field", each FORMAT field is provided as a separate column with nested
        sample name fields.
    samples : ``"*"``, list[str], or None, optional [default: ``None``]
        Samples to include in the genotype output. ``"*"`` includes all samples
        declared in the header. Pass a list to select specific samples. ``None``
        omits all sample genotype data.
    samples_nested : bool, optional [default: False]
        Whether to nest sample data under a single structured column.
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
        genotype_fields=genotype_fields,
        genotype_by=genotype_by,
        samples=samples,
        samples_nested=samples_nested,
        coords=coords,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )


def from_bcf(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    compression: Literal["bgzf", None] = "bgzf",
    *,
    fields: Literal["*"] | list[str] | None = "*",
    info_fields: Literal["*"] | list[str] | None = "*",
    genotype_fields: Literal["*"] | list[str] | None = "*",
    genotype_by: Literal["sample", "field"] = "sample",
    samples: Literal["*"] | list[str] | None = None,
    samples_nested: bool = False,
    coords: Literal["01", "11"] = "11",
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> BcfFile:
    """
    Create a BCF file data source.

    .. versionchanged:: 0.7.0
        The ``samples`` parameter now defaults to omitting sample genotype
        data (``None``) instead of including all samples (``"*"``). To include
        samples, pass a value to the ``samples`` parameter or use the
        :meth:`~oxbow.core.BcfFile.with_samples()` method on the returned data
        source.

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
    fields : ``"*"``, list[str], or None, optional [default: ``"*"``]
        Fixed fields to project. ``"*"`` includes all standard fields. Pass a
        list to select specific fields. ``None`` omits all fixed fields.
    info_fields : ``"*"``, list[str], or None, optional [default: ``"*"``]
        INFO fields to project, nested under an ``"info"`` column. ``"*"``
        includes all INFO fields declared in the header. Pass a list to select
        specific fields. ``None`` omits the info column entirely.
    genotype_fields : ``"*"``, list[str], or None, optional [default: ``"*"``]
        Genotype (aka "FORMAT") fields to project for each sample. ``"*"``
        includes all FORMAT fields declared in the header. Pass a list to select
        specific fields. ``None`` omits the genotype fields.
    genotype_by : Literal["sample", "field"], optional [default: "sample"]
        Determines how genotype-specific data is organized. If "sample", each
        sample is provided as a separate column with nested FORMAT fields. If
        "field", each FORMAT field is provided as a separate column with nested
        sample name fields.
    samples : ``"*"``, list[str], or None, optional [default: ``None``]
        Samples to include in the genotype output. ``"*"`` includes all samples
        declared in the header. Pass a list to select specific samples. ``None``
        omits all sample genotype data.
    samples_nested : bool, optional [default: False]
        Whether to nest sample data under a single structured column.
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
        genotype_fields=genotype_fields,
        genotype_by=genotype_by,
        samples=samples,
        samples_nested=samples_nested,
        coords=coords,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )
