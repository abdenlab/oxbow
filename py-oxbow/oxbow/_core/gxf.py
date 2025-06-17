"""
DataSource classes for GTF/GFF3 formats.
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
from oxbow.oxbow import PyGffScanner, PyGtfScanner


class GxfFile(DataSource):
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
                if "attributes" not in columns:
                    scan_kwargs["attribute_defs"] = []
                elif scan_kwargs.get("attribute_defs") == []:
                    raise ValueError(
                        "Cannot select `attributes` column if no attribute "
                        "definitions are provided."
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
        attribute_defs: list[tuple[str, str]] | None = None,
        attribute_scan_rows: int = 1024,
        regions: str | list[str] | None = None,
        index: str | Callable[[], IO[bytes] | str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        super().__init__(source, index, batch_size)

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._scanner_kwargs = dict(compressed=compressed)
        if attribute_defs is None:
            attribute_defs = self.scanner().attribute_defs(attribute_scan_rows)
        self._schema_kwargs = dict(fields=fields, attribute_defs=attribute_defs)

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
    def attribute_defs(self) -> list[tuple[str, str]]:
        """List of definitions for interpreting attribute records."""
        return self._schema_kwargs["attribute_defs"]


class GtfFile(GxfFile):
    _scanner_type = PyGtfScanner


class GffFile(GxfFile):
    _scanner_type = PyGffScanner


def from_gtf(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    compression: Literal["infer", "bgzf", "gzip", None] = "infer",
    *,
    fields: list[str] | None = None,
    attribute_defs: list[tuple[str, str]] | None = None,
    attribute_scan_rows: int = 1024,
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> GtfFile:
    """
    Create a GTF file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the GTF file, or a callable that opens the file
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
    attribute_defs : list[tuple[str, str]], optional [default: None]
        Definitions for variable attribute fields to project. These will be
        nested in an "attributes" column. If None, attribute definitions are
        discovered by scanning records in the file, which is controlled by the
        ``attribute_scan_rows`` parameter. To omit attributes entirely,
        set ``attribute_defs=[]``.
    attribute_scan_rows : int, optional [default: 1024]
        Number of rows to scan for attribute definitions.
    regions : str | list[str], optional
        One or more genomic regions to query. Only applicable if an associated
        index file is available.
    index : str, pathlib.Path, or Callable, optional
        An optional index file associated with the GTF file. If ``source`` is a
        URI or path, is BGZF-compressed, and the index file shares the same
        name with a ".tbi" or ".csi" extension, the index file is automatically
        detected.
    batch_size : int, optional [default: 131072]
        The number of records to read in each batch.

    Returns
    -------
    GtfFile
        A data source object representing the GTF file.

    See also
    --------
    from_gff : Create a GFF file data source.
    from_bed : Create a BED file data source.
    from_bigbed : Create a BigBed file data source.
    from_bigwig : Create a BigWig file data source.
    """
    source, index, bgzf_compressed = prepare_source_and_index(
        source, index, compression
    )
    return GtfFile(
        source=source,
        compressed=bgzf_compressed,
        fields=fields,
        attribute_defs=attribute_defs,
        attribute_scan_rows=attribute_scan_rows,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )


def from_gff(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    compression: Literal["infer", "bgzf", "gzip", None] = "infer",
    *,
    fields: list[str] | None = None,
    attribute_defs: list[tuple[str, str]] | None = None,
    attribute_scan_rows: int = 1024,
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> GffFile:
    """
    Create a GFF3 file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the GFF file, or a callable that opens the file
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
    attribute_defs : list[tuple[str, str]], optional [default: None]
        Definitions for variable attribute fields to project. These will be
        nested in an "attributes" column. If None, attribute definitions are
        discovered by scanning records in the file, which is controlled by the
        ``attribute_scan_rows`` parameter. To omit attributes entirely,
        set ``attribute_defs=[]``.
    attribute_scan_rows : int, optional [default: 1024]
        Number of rows to scan for attribute definitions.
    regions : str | list[str], optional
        One or more genomic regions to query. Only applicable if an associated
        index file is available.
    index : str, pathlib.Path, or Callable, optional
        An optional index file associated with the GTF file. If ``source`` is a
        URI or path, is BGZF-compressed, and the index file shares the same
        name with a ".tbi" or ".csi" extension, the index file is automatically
        detected.
    batch_size : int, optional [default: 131072]
        The number of records to read in each batch.

    Returns
    -------
    GffFile
        A data source object representing the GFF file.

    See also
    --------
    from_gtf : Create a GTF file data source.
    from_bed : Create a BED file data source.
    from_bigbed : Create a BigBed file data source.
    from_bigwig : Create a BigWig file data source.
    """
    source, index, bgzf_compressed = prepare_source_and_index(
        source, index, compression
    )
    return GffFile(
        source=source,
        fields=fields,
        index=index,
        compressed=bgzf_compressed,
        regions=regions,
        attribute_defs=attribute_defs,
        attribute_scan_rows=attribute_scan_rows,
        batch_size=batch_size,
    )
