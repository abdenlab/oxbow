"""
DataSource classes for htslib alignment formats.
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
from oxbow.oxbow import PyBamScanner, PySamScanner


class AlignmentFile(DataSource):
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
                if "tags" not in columns:
                    scan_kwargs["tag_defs"] = []
                elif scan_kwargs.get("tag_defs") == []:
                    raise ValueError(
                        "Cannot select `tags` column if no tag definitions are "
                        "provided."
                    )

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
                if region == "*":
                    yield self._batchreader_builder(
                        scanner.scan_unmapped, scanner.field_names()
                    )
                else:
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
        fields: list[str] | None = None,
        tag_defs: list[tuple[str, str]] | None = None,
        tag_scan_rows: int = 1024,
        regions: str | list[str] | None = None,
        index: str | Callable[[], IO[bytes] | str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        super().__init__(source, index, batch_size)

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._scanner_kwargs = dict(compressed=compressed)
        if tag_defs is None:
            tag_defs = self.scanner().tag_defs(tag_scan_rows)
        self._schema_kwargs = dict(fields=fields, tag_defs=tag_defs)

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
    def tag_defs(self) -> list[tuple[str, str]]:
        """List of definitions for interpreting tag records."""
        return self._schema_kwargs["tag_defs"]


class SamFile(AlignmentFile):
    _scanner_type = PySamScanner


class BamFile(AlignmentFile):
    _scanner_type = PyBamScanner


def from_sam(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    compression: Literal["infer", "bgzf", "gzip", None] = "infer",
    *,
    fields: list[str] | None = None,
    tag_defs: list[tuple[str, str]] = None,
    tag_scan_rows: int = 1024,
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> SamFile:
    """
    Create a SAM file data source.

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
    fields : list[str], optional
        Specific fixed fields to project. By default, all fixed fields are
        included.
    tag_defs : list[tuple[str, str]], optional [default: None]
        Definitions for variable tag fields to project. These will be nested in
        a "tags" column. If None, tag definitions are discovered by scanning
        records in the file, which is controlled by the ``tag_scan_rows``
        parameter. To omit tags entirely, set ``tag_defs=[]``.
    tag_scan_rows : int, optional [default: 1024]
        Number of rows to scan for tag definitions.
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
    """
    source, index, bgzf_compressed = prepare_source_and_index(
        source, index, compression
    )
    return SamFile(
        source=source,
        compressed=bgzf_compressed,
        fields=fields,
        tag_defs=tag_defs,
        tag_scan_rows=tag_scan_rows,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )


def from_bam(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    compression: Literal["bgzf", None] = "bgzf",
    *,
    fields: list[str] | None = None,
    tag_defs: list[tuple[str, str]] = None,
    tag_scan_rows: int = 1024,
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> BamFile:
    """
    Create a BAM file data source.

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
    fields : list[str], optional
        Specific fixed fields to project. By default, all fixed fields are
        included.
    tag_defs : list[tuple[str, str]], optional [default: None]
        Definitions for variable tag fields to project. These will be nested in
        a "tags" column. If None, tag definitions are discovered by scanning
        records in the file, which is controlled by the ``tag_scan_rows``
        parameter. To omit tags entirely, set ``tag_defs=[]``.
    tag_scan_rows : int, optional [default: 1024]
        Number of rows to scan for tag definitions.
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
    BamFile
        A data source object representing the BAM file.

    Notes
    -----
    Binary Alignment Map (BAM) is a binary representation of SAM files.

    See also
    --------
    from_sam : Create a SAM file data source.
    """
    source, index, bgzf_compressed = prepare_source_and_index(
        source, index, compression
    )
    return BamFile(
        source=source,
        compressed=bgzf_compressed,
        fields=fields,
        tag_defs=tag_defs,
        tag_scan_rows=tag_scan_rows,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )
