"""
DataSource classes for the BED family of formats.
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
from oxbow.oxbow import PyBedScanner


class BedFile(DataSource):
    _scanner_type = PyBedScanner

    def _batchreader_builder(
        self,
        scan_fn: Callable,
        field_names: list[str],
        region: str | None = None,
    ) -> Callable[[list[str] | None, int], pa.RecordBatchReader]:
        def builder(columns, batch_size):
            scan_kwargs = self._schema_kwargs.copy()

            # TODO: Subsetting BED fields is not supported yet.
            # field_names = [name for name in self._scanner.field_names()]
            # if columns is not None:
            #     scan_kwargs["fields"] = [col for col in columns if col in field_names]

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
        bed_schema: str = "bed3+",
        compressed: bool = False,
        *,
        fields: list[str] | None = None,
        regions: str | list[str] | None = None,
        index: str | Callable[[], IO[bytes] | str] | None = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ):
        super().__init__(source, index, batch_size)

        if isinstance(regions, str):
            regions = [regions]
        self._regions = regions

        self._scanner_kwargs = dict(bed_schema=bed_schema, compressed=compressed)
        self._schema_kwargs = dict(fields=fields)

    def regions(self, regions: str | list[str]) -> Self:
        return type(self)(
            self._src,
            regions=regions,
            index=self._index_src,
            batch_size=self._batch_size,
            **self._scanner_kwargs,
            **self._schema_kwargs,
        )


def from_bed(
    source: str | pathlib.Path | Callable[[], IO[bytes] | str],
    bed_schema: str = "bed3+",
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
    bed_schema : str, optional [default: "bed3+"]
        Schema for intepreting the BED file. The default is "bed3+", which
        includes the first three standard fields (chrom, start, end) and
        any additional data is lumped into a single "rest" column.
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
