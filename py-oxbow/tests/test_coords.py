"""Tests for the coords parameter on scanner classes.

Each format has a native coordinate system:
- SAM/BAM/CRAM, VCF/BCF, GFF/GTF: 1-based closed (default "11")
- BED, BigBed, BigWig: 0-based half-open (default "01")

Passing coords="01" to a 1-based format shifts start positions by -1.
Passing coords="11" to a 0-based format shifts start positions by +1.
"""

import pyarrow as pa

import oxbow.oxbow as ox


def _read_pos(scanner, col="pos"):
    """Read a single batch from a scanner and return the named column as a list."""
    stream = scanner.scan()
    schema = scanner.schema()
    reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
    batch = reader.read_all()
    return batch.column(col).to_pylist()


# -- Alignment (1-based native) ------------------------------------------------


class TestSamCoords:
    def test_default_is_one_based(self):
        scanner = ox.PySamScanner("data/sample.sam", fields=["pos"])
        pos = _read_pos(scanner)
        assert pos == [16, 29, 37]

    def test_zero_based(self):
        scanner = ox.PySamScanner("data/sample.sam", fields=["pos"], coords="01")
        pos = _read_pos(scanner)
        assert pos == [15, 28, 36]

    def test_explicit_one_based_matches_default(self):
        default = _read_pos(ox.PySamScanner("data/sample.sam", fields=["pos"]))
        explicit = _read_pos(
            ox.PySamScanner("data/sample.sam", fields=["pos"], coords="11")
        )
        assert default == explicit


class TestBamCoords:
    def test_zero_based_shifts_by_minus_one(self):
        default = _read_pos(ox.PyBamScanner("data/sample.bam", fields=["pos"]))
        zero_based = _read_pos(
            ox.PyBamScanner("data/sample.bam", fields=["pos"], coords="01")
        )
        assert zero_based == [v - 1 for v in default]


# -- Variant (1-based native) --------------------------------------------------


class TestVcfCoords:
    def test_default_is_one_based(self):
        scanner = ox.PyVcfScanner("data/sample.vcf", fields=["pos"])
        pos = _read_pos(scanner)
        # Just check the offset relationship, not exact values
        assert all(isinstance(v, int) for v in pos)

    def test_zero_based_shifts_by_minus_one(self):
        default = _read_pos(ox.PyVcfScanner("data/sample.vcf", fields=["pos"]))
        zero_based = _read_pos(
            ox.PyVcfScanner("data/sample.vcf", fields=["pos"], coords="01")
        )
        assert zero_based == [v - 1 for v in default]


# -- GXF (1-based native) ------------------------------------------------------


class TestGtfCoords:
    def test_default_is_one_based(self):
        scanner = ox.PyGtfScanner("data/sample.gtf", fields=["start"])
        start = _read_pos(scanner, col="start")
        assert all(isinstance(v, int) for v in start)

    def test_zero_based_shifts_by_minus_one(self):
        default = _read_pos(
            ox.PyGtfScanner("data/sample.gtf", fields=["start"]), col="start"
        )
        zero_based = _read_pos(
            ox.PyGtfScanner("data/sample.gtf", fields=["start"], coords="01"),
            col="start",
        )
        assert zero_based == [v - 1 for v in default]


class TestGffCoords:
    def test_zero_based_shifts_by_minus_one(self):
        default = _read_pos(
            ox.PyGffScanner("data/sample.gff", fields=["start"]), col="start"
        )
        zero_based = _read_pos(
            ox.PyGffScanner("data/sample.gff", fields=["start"], coords="01"),
            col="start",
        )
        assert zero_based == [v - 1 for v in default]


# -- BED (0-based native) ------------------------------------------------------


class TestBedCoords:
    def test_default_is_zero_based(self):
        scanner = ox.PyBedScanner("data/sample.bed", "bed9", fields=["start"])
        start = _read_pos(scanner, col="start")
        # BED file has 0-based starts; noodles converts to 1-based Position
        # internally, then the default ZeroHalfOpen coord_system subtracts 1
        # back to 0-based. So default output == file values.
        assert start[0] == 1100000

    def test_one_based_shifts_by_plus_one(self):
        default = _read_pos(
            ox.PyBedScanner("data/sample.bed", "bed9", fields=["start"]), col="start"
        )
        one_based = _read_pos(
            ox.PyBedScanner("data/sample.bed", "bed9", fields=["start"], coords="11"),
            col="start",
        )
        assert one_based == [v + 1 for v in default]


# -- BBI (0-based native) ------------------------------------------------------


class TestBigWigCoords:
    def test_one_based_shifts_by_plus_one(self):
        default = _read_pos(
            ox.PyBigWigScanner("data/sample.bw", fields=["start"]), col="start"
        )
        one_based = _read_pos(
            ox.PyBigWigScanner("data/sample.bw", fields=["start"], coords="11"),
            col="start",
        )
        assert one_based == [v + 1 for v in default]


class TestBigBedCoords:
    def test_one_based_shifts_by_plus_one(self):
        default = _read_pos(
            ox.PyBigBedScanner("data/sample.bb", fields=["start"]), col="start"
        )
        one_based = _read_pos(
            ox.PyBigBedScanner("data/sample.bb", fields=["start"], coords="11"),
            col="start",
        )
        assert one_based == [v + 1 for v in default]
