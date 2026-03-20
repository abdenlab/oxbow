"""
DataSource classes for htslib alignment formats.
"""

from __future__ import annotations

import pathlib
from typing import IO, Callable, Literal

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

from oxbow._core.base import DEFAULT_BATCH_SIZE, DataSource, prepare_source_and_index
from oxbow.oxbow import PyBamScanner, PySamScanner, PyCramScanner


class AlignmentFile(DataSource):
    def __init__(
        self,
        source: str | Callable[[], IO[bytes] | str],
        compressed: bool = False,
        *,
        fields: Literal["*"] | list[str] | None = "*",
        tag_defs: list[tuple[str, str]] | None = None,
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
            compressed=compressed, fields=fields, tag_defs=tag_defs, coords=coords
        )

    def _tag_discovery_kwargs(self) -> dict:
        return dict(compressed=self._scanner_kwargs["compressed"])

    def _scan_query(self, scanner, region, columns, batch_size):
        if region == "*":
            return scanner.scan_unmapped(
                index=self._index, columns=columns, batch_size=batch_size
            )
        else:
            return scanner.scan_query(
                region=region, index=self._index, columns=columns, batch_size=batch_size
            )

    def with_tags(
        self, tag_defs: list[tuple[str, str]] | None = None, *, scan_rows: int = 1024
    ) -> Self:
        """
        Return a new data source with the specified tag definitions.

        Parameters
        ----------
        tag_defs : list[tuple[str, str]] or None, optional [default: None]
            Definitions for tags to project. These will be nested in a "tags"
            column. If None (default), tag definitions are discovered by
            scanning records in the file, which is controlled by the
            ``scan_rows`` parameter.
        scan_rows : int, optional [default: 1024]
            Number of rows to scan for tag discovery if tag_defs is None. Set
            to -1 to scan the entire file, which may be slow for large files.

        Returns
        -------
        Self
            A new data source with the specified tag definitions.

        Notes
        -----
        Tag definitions take the form of a list of (tag_name, tag_type) tuples,
        where tag_name is a 2-character string and tag_type is a
        single-character type code as defined in the SAM specification.

        Type codes:

            - A: Printable character
            - i: Signed integer
            - f: Floating point number
            - Z: String
            - H: Hex string
            - B: Array (comma-separated values with type code prefix, e.g., "i,1,2,3")
        """
        if tag_defs is None:
            scan_rows = scan_rows if scan_rows >= 0 else None
            discovered = self._scanner_type(
                self._source, **self._tag_discovery_kwargs()
            ).tag_defs(scan_rows)
            self._scanner_kwargs["tag_defs"] = discovered or []

        return type(self)(
            self._src,
            regions=self._regions,
            index=self._index_src,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
        )

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._src,
            regions=regions,
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
    def tag_defs(self) -> list[tuple[str, str]]:
        """List of definitions for interpreting tags."""
        return self._scanner_kwargs["tag_defs"]


class SamFile(AlignmentFile):
    _scanner_type = PySamScanner


class BamFile(AlignmentFile):
    _scanner_type = PyBamScanner


class CramFile(AlignmentFile):
    _scanner_type = PyCramScanner

    def __init__(
        self,
        source: str | Callable[[], IO[bytes] | str],
        compressed: bool = False,
        *,
        fields: Literal["*"] | list[str] | None = "*",
        tag_defs: list[tuple[str, str]] | None = None,
        coords: Literal["01", "11"] = "11",
        regions: str | list[str] | None = None,
        index: str | Callable[[], IO[bytes] | str] | None = None,
        reference: str | Callable[[], IO[bytes] | str] | None = None,
        reference_index: str | Callable[[], IO[bytes] | str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        self._reference = reference
        self._reference_index = reference_index
        super().__init__(
            source,
            compressed=compressed,
            fields=fields,
            tag_defs=tag_defs,
            coords=coords,
            regions=regions,
            index=index,
            batch_size=batch_size,
        )
        self._scanner_kwargs["reference"] = reference
        self._scanner_kwargs["reference_index"] = reference_index

    def _tag_discovery_kwargs(self) -> dict:
        return dict(
            compressed=self._scanner_kwargs.get("compressed", False),
            reference=self._reference,
            reference_index=self._reference_index,
        )


def from_sam(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    compression: Literal["infer", "bgzf", "gzip", None] = "infer",
    *,
    fields: Literal["*"] | list[str] | None = "*",
    tag_defs: list[tuple[str, str]] | None = None,
    coords: Literal["01", "11"] = "11",
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> SamFile:
    """
    Create a SAM file data source.

    .. versionchanged:: 0.7.0
        The ``tag_scan_rows`` parameter was removed and tag definitions are no
        longer discovered by default. The ``tag_defs`` parameter now defaults
        to omitting tag definitions (``None``). To perform tag discovery,
        use the :meth:`~oxbow.core.SamFile.with_tags()` method on the returned
        data source, which accepts a ``scan_rows`` parameter to control how
        many records are scanned.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the SAM file, or a callable that opens the file
        as a file-like object.
    compression : Literal["infer", "bgzf", "gzip", None], default: "infer"
        Compression of the source bytestream. If "infer" and ``source`` is a
        URI or path, the file's compression is guessed based on the extension,
        where ".gz" or ".bgz" is interpreted as BGZF. Pass "gzip" to decode
        regular GZIP. If None, the source bytestream is assumed to be
        uncompressed. For more customized decoding, provide a callable
        ``source`` instead.
    fields : list[str] or "*", optional [default: "*"]
        Standard SAM fields to include. By default, all standard fields are
        included.
    tag_defs : list[tuple[str, str]], optional [default: None]
        Definitions for tags to project. These will be nested in a "tags"
        column. If None, tag definitions are omitted. To discover tag
        definitions, use the ``with_tags()`` method on the returned data source.
    regions : str | list[str], optional
        One or more genomic regions to query. Only applicable if an associated
        index file is available.
    index : str, pathlib.Path, or Callable, optional
        An optional index file associated with the SAM file. If ``source`` is a
        URI or path, is BGZF-compressed, and the index file shares the same
        name with a ".tbi" or ".csi" extension, the index file is automatically
        detected.
    batch_size : int, optional [default: 131072]
        The number of records to read in each batch.

    Returns
    -------
    SamFile
        A data source object representing the SAM file.

    Notes
    -----
    Sequence Alignment Map (SAM) is a widely used text-based format for
    storing biological sequences aligned to a reference sequence.

    See also
    --------
    from_bam : Create a BAM file data source.
    from_cram : Create a CRAM file data source.
    """
    source, index, bgzf_compressed = prepare_source_and_index(
        source, index, compression
    )
    return SamFile(
        source=source,
        compressed=bgzf_compressed,
        fields=fields,
        tag_defs=tag_defs,
        coords=coords,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )


def from_bam(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    compression: Literal["bgzf", None] = "bgzf",
    *,
    fields: Literal["*"] | list[str] | None = "*",
    tag_defs: list[tuple[str, str]] | None = None,
    coords: Literal["01", "11"] = "11",
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> BamFile:
    """
    Create a BAM file data source.

    .. versionchanged:: 0.7.0
        The ``tag_scan_rows`` parameter was removed and tag definitions are no
        longer discovered by default. The ``tag_defs`` parameter now defaults
        to omitting tag definitions (``None``). To perform tag discovery,
        use the :meth:`~oxbow.core.BamFile.with_tags()` method on the returned
        data source, which accepts a ``scan_rows`` parameter to control how
        many records are scanned.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the BAM file, or a callable that opens the file
        as a file-like object.
    compression : Literal["bgzf", None], default: "bgzf"
        Compression of the source bytestream. By default, BAM sources are
        assumed to be BGZF-compressed. If None, the source is assumed to be
        uncompressed. For more custom decoding, provide a callable ``source``
        instead.
    fields : list[str] or "*", optional [default: "*"]
        Standard SAM fields to include. By default, all standard fields are
        included.
    tag_defs : list[tuple[str, str]], optional [default: None]
        Definitions for tags to project. These will be nested in a "tags"
        column. If None, tag definitions are omitted. To discover tag
        definitions, use the ``with_tags()`` method on the returned data source.
    regions : str | list[str], optional
        One or more genomic regions to query. Only applicable if an associated
        index file is available.
    index : str, pathlib.Path, or Callable, optional
        An optional index file associated with the BAM file. If ``source`` is a
        URI or path, is BGZF-compressed, and the index file shares the same
        name with a ".tbi" or ".csi" extension, the index file is automatically
        detected.
    batch_size : int, optional [default: 131072]
        The number of records to read in each batch.

    Returns
    -------
    BamFile
        A data source object representing the BAM file.

    Notes
    -----
    Binary Alignment Map (BAM) is a binary representation of SAM files.

    See also
    --------
    from_sam : Create a SAM file data source.
    from_cram : Create a CRAM file data source.
    """
    source, index, bgzf_compressed = prepare_source_and_index(
        source, index, compression
    )
    return BamFile(
        source=source,
        compressed=bgzf_compressed,
        fields=fields,
        tag_defs=tag_defs,
        coords=coords,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )


def from_cram(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    *,
    fields: Literal["*"] | list[str] | None = "*",
    tag_defs: list[tuple[str, str]] | None = None,
    coords: Literal["01", "11"] = "11",
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    reference: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    reference_index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> CramFile:
    """
    Create a CRAM file data source.

    .. versionchanged:: 0.7.0
        The ``tag_scan_rows`` parameter was removed and tag definitions are no
        longer discovered by default. The ``tag_defs`` parameter now defaults
        to omitting tag definitions (``None``). To perform tag discovery,
        use the :meth:`~oxbow.core.CramFile.with_tags()` method on the returned
        data source, which accepts a ``scan_rows`` parameter to control how
        many records are scanned.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the CRAM file, or a callable that opens the file
        as a file-like object.
    fields : list[str] or "*", optional [default: "*"]
        Standard SAM fields to include. By default, all standard fields are
        included.
    tag_defs : list[tuple[str, str]], optional [default: None]
        Definitions for tags to project. These will be nested in a "tags"
        column. If None, tag definitions are omitted. To discover tag
        definitions, use the ``with_tags()`` method on the returned data source.
    regions : str | list[str], optional
        One or more genomic regions to query. Only applicable if an associated
        index file is available.
    index : str, pathlib.Path, or Callable, optional
        An optional index file associated with the CRAM file. If ``source`` is
        a URI or path and the index file shares the same name with a ".crai"
        the index file is automatically detected.
    reference : str, pathlib.Path, or Callable, optional
        The URI or path to the FASTA reference file used for CRAM encoding, or
        a callable that opens the file as a file-like object. Required if the
        CRAM file does not contain an embedded reference.
    reference_index : str, pathlib.Path, or Callable, optional
        The URI or path to the FASTA reference index file, or a callable that
        opens the file as a file-like object. If ``reference`` is provided as a
        URI or path and the index file shares the same name with a ".fai"
        extension, the index file is automatically detected.
    batch_size : int, optional [default: 131072]
        The number of records to read in each batch.

    Returns
    -------
    CramFile
        A data source object representing the CRAM file.

    Notes
    -----
    CRAM is a compressed binary format for storing sequence alignments.

    See also
    --------
    from_sam : Create a SAM file data source.
    from_bam : Create a BAM file data source.
    """
    source, index, _ = prepare_source_and_index(source, index, False)
    return CramFile(
        source=source,
        fields=fields,
        tag_defs=tag_defs,
        coords=coords,
        regions=regions,
        index=index,
        reference=reference,
        reference_index=reference_index,
        batch_size=batch_size,
    )
