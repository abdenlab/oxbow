"""
DataSource classes for the BED family of formats.
"""

from __future__ import annotations

import pathlib
from typing import IO, TYPE_CHECKING, Callable, Literal, Union

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

from oxbow._core.base import DEFAULT_BATCH_SIZE, DataSource, prepare_source_and_index
from oxbow.oxbow import PyBedScanner

if TYPE_CHECKING:
    CustomFieldDefs = Union[list[tuple[str, str]], dict[str, str]]
    BedSchemaLike = Union[str, tuple[str, CustomFieldDefs]]


class BedFile(DataSource):
    _scanner_type = PyBedScanner

    def __init__(
        self,
        source: str | Callable[[], IO[bytes] | str],
        bed_schema: BedSchemaLike = "bed3+",
        compressed: bool = False,
        *,
        fields: Literal["*"] | list[str] | None = "*",
        regions: str | list[str] | None = None,
        index: str | Callable[[], IO[bytes] | str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        super().__init__(source, index, batch_size)

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._scanner_kwargs = dict(
            bed_schema=bed_schema,
            compressed=compressed,
            fields=fields,
        )

    def _scan_query(self, scanner, region, columns, batch_size):
        return scanner.scan_query(
            region=region,
            index=self._index,
            columns=columns,
            batch_size=batch_size,
        )

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._src,
            regions=regions,
            index=self._index_src,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
        )


def from_bed(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    bed_schema: BedSchemaLike = "bed3+",
    compression: Literal["infer", "bgzf", "gzip", None] = "infer",
    *,
    fields: list[str] | None = None,
    regions: str | list[str] | None = None,
    index: str | pathlib.Path | Callable[[], IO[bytes] | str] | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> BedFile:
    """
    Create a BED file data source.

    Parameters
    ----------
    source : str, pathlib.Path, or Callable
        The URI or path to the BED file, or a callable that opens the file
        as a file-like object.
    bed_schema : str or tuple[str, list | dict], optional
        Schema for interpreting the BED file. Can be a specifier string
        (e.g., "bed6", "bed3+", "bedgraph") or a tuple of
        ``(base_specifier, custom_defs)`` where ``custom_defs`` is a list
        of ``(name, type)`` tuples or a dict mapping field names to type
        strings. The default is "bed3+", which includes the first three
        standard fields and lumps any additional data into a "rest" column.
    compression : Literal["infer", "bgzf", "gzip", None], default: "infer"
        Compression of the source bytestream. If "infer" and ``source`` is a
        URI or path, the file's compression is guessed based on the extension,
        where ".gz" or ".bgz" is interpreted as BGZF. Pass "gzip" to decode
        regular GZIP. If None, the source bytestream is assumed to be
        uncompressed. For more customized decoding, provide a callable
        ``source`` instead.
    fields : list[str], optional
        Specific fields to project as columns. By default, all available fields
        are included.
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
    BedFile
        A data source object representing the BED file.

    See also
    --------
    from_gtf : Create a GTF file data source.
    from_gff : Create a GFF file data source.
    from_bigbed : Create a BigBed file data source.
    from_bigwig : Create a BigWig file data source.
    """
    source, index, bgzf_compressed = prepare_source_and_index(
        source, index, compression
    )
    return BedFile(
        source=source,
        bed_schema=bed_schema,
        compressed=bgzf_compressed,
        fields=fields,
        regions=regions,
        index=index,
        batch_size=batch_size,
    )
