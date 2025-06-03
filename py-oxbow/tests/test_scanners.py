import pickle

import pyarrow as pa
import pytest
import pytest_manifest

import oxbow.oxbow as ox
from tests.utils import Input


class TestPySamScanner:
    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("qname", "rname", "mapq")),
        ],
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(pickle.loads(pickle.dumps(ox.PySamScanner("data/sample.sam"))))
        )
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("qname", "rname", "foo"))
        scanner = pickle.loads(
            pickle.dumps(pickle.loads(pickle.dumps(ox.PySamScanner("data/sample.sam"))))
        )
        error = None
        try:
            scanner.scan(*input.args, **input.kwargs)
        except ValueError as e:
            error = str(e)
        finally:
            assert manifest == error

    def test_pickle(self):
        scanner = pickle.loads(pickle.dumps(ox.PySamScanner("data/sample.sam")))
        assert isinstance(scanner, ox.PySamScanner)


class TestPyBamScanner:
    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("qname", "rname", "mapq")),
        ],
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(pickle.loads(pickle.dumps(ox.PyBamScanner("data/sample.bam"))))
        )
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("qname", "rname", "foo"))
        scanner = pickle.loads(
            pickle.dumps(pickle.loads(pickle.dumps(ox.PyBamScanner("data/sample.bam"))))
        )
        error = None
        try:
            scanner.scan(*input.args, **input.kwargs)
        except ValueError as e:
            error = str(e)
        finally:
            assert manifest == error

    def test_pickle(self):
        scanner = pickle.loads(pickle.dumps(ox.PyBamScanner("data/sample.bam")))
        assert isinstance(scanner, ox.PyBamScanner)


class TestPyBcfScanner:
    def test_chrom_names(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyBcfScanner("data/sample.bcf")))
        assert manifest == scanner.chrom_names()

    def test_chrom_sizes(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyBcfScanner("data/sample.bcf")))
        assert manifest == scanner.chrom_sizes()

    def test_info_field_defs(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyBcfScanner("data/sample.bcf")))
        assert manifest == scanner.info_field_defs()

    def test_genotype_field_defs(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyBcfScanner("data/sample.bcf")))
        assert manifest == scanner.genotype_field_defs()

    def test_sample_names(self):
        scanner = pickle.loads(pickle.dumps(ox.PyBcfScanner("data/sample.bcf")))
        assert 1233 == len(scanner.sample_names())

    @pytest.mark.parametrize(
        "input",
        Input.permute(batch_size=[2], samples=[["HG00096", "HG00101", "HG00103"]]),
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyBcfScanner("data/sample.bcf")))
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("name", "sequence", "foo"))
        scanner = pickle.loads(pickle.dumps(ox.PyBcfScanner("data/sample.bcf")))
        error = None
        try:
            scanner.scan(*input.args, **input.kwargs)
        except ValueError as e:
            error = str(e)
        finally:
            assert manifest == error

    @pytest.mark.parametrize(
        "input",
        Input.permute(
            region=["Y"],
            index=["data/sample.bcf.csi"],
            batch_size=[2],
            samples=[["HG00096", "HG00101", "HG00103"]],
        ),
    )
    def test_scan_query(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyBcfScanner("data/sample.bcf")))
        schema = scanner.schema()
        stream = scanner.scan_query(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_pickle(self):
        scanner = pickle.loads(pickle.dumps(ox.PyBcfScanner("data/sample.bcf")))
        assert isinstance(scanner, ox.PyBcfScanner)


class TestPyVcfScanner:
    def test_chrom_names(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf", compressed=False))
        )
        assert manifest == scanner.chrom_names()

    def test_chrom_names_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf.gz", compressed=True))
        )
        assert manifest == scanner.chrom_names()

    def test_chrom_sizes(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf", compressed=False))
        )
        assert manifest == scanner.chrom_sizes()

    def test_chrom_sizes_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf.gz", compressed=True))
        )
        assert manifest == scanner.chrom_sizes()

    def test_info_field_names(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf", compressed=False))
        )
        assert manifest == scanner.info_field_names()

    def test_info_field_names_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf.gz", compressed=True))
        )
        assert manifest == scanner.info_field_names()

    def test_info_field_defs(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf", compressed=False))
        )
        assert manifest == scanner.info_field_defs()

    def test_info_field_defs_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf.gz", compressed=True))
        )
        assert manifest == scanner.info_field_defs()

    def test_genotype_field_names(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf", compressed=False))
        )
        assert manifest == scanner.genotype_field_names()

    def test_genotype_field_names_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf.gz", compressed=True))
        )
        assert manifest == scanner.genotype_field_names()

    def test_genotype_field_defs(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf", compressed=False))
        )
        assert manifest == scanner.genotype_field_defs()

    def test_genotype_field_defs_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf.gz", compressed=True))
        )
        assert manifest == scanner.genotype_field_defs()

    def test_sample_names(self):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf", compressed=False))
        )
        assert 3 == len(scanner.sample_names())

    def test_sample_names_compressed(self):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf.gz", compressed=True))
        )
        assert 3 == len(scanner.sample_names())

    @pytest.mark.parametrize(
        "input",
        Input.permute(batch_size=[2], samples=[["NA12878i", "NA12891", "NA12892"]]),
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf", compressed=False))
        )
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    @pytest.mark.parametrize(
        "input",
        Input.permute(batch_size=[2], samples=[["NA12878i", "NA12891", "NA12892"]]),
    )
    def test_scan_compressed(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf.gz", compressed=True))
        )
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("name", "sequence", "foo"))
        scanner = pickle.loads(pickle.dumps(ox.PyVcfScanner("data/sample.vcf")))
        error = None
        try:
            scanner.scan(*input.args, **input.kwargs)
        except ValueError as e:
            error = str(e)
        finally:
            assert manifest == error

    @pytest.mark.parametrize(
        "input",
        Input.permute(
            region=["X"],
            index=["data/sample.vcf.gz.csi", "data/sample.vcf.gz.tbi"],
            batch_size=[2],
            samples=[["NA12878i", "NA12891", "NA12892"]],
        ),
    )
    def test_scan_query(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf.gz", compressed=True))
        )
        schema = scanner.schema()
        stream = scanner.scan_query(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_pickle(self):
        scanner = pickle.loads(
            pickle.dumps(ox.PyVcfScanner("data/sample.vcf", compressed=False))
        )
        assert isinstance(scanner, ox.PyVcfScanner)


class TestPyFastaScanner:
    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("name", "sequence")),
        ],
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyFastaScanner("data/sample.fasta")))
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    @pytest.mark.parametrize(
        "input",
        [
            Input(regions=["seq1:10-20", "seq10"], index="data/sample.fasta.fai"),
            Input(
                regions=["seq1:10-20", "seq10", "seq2:1-30"],
                index="data/sample.fasta.fai",
            ),
        ],
    )
    def test_scan_query(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyFastaScanner("data/sample.fasta")))
        schema = scanner.schema()
        stream = scanner.scan_query(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("name", "sequence", "foo"))
        scanner = pickle.loads(pickle.dumps(ox.PyFastaScanner("data/sample.fasta")))
        error = None
        try:
            scanner.scan(*input.args, **input.kwargs)
        except ValueError as e:
            error = str(e)
        finally:
            assert manifest == error

    def test_pickle(self):
        scanner = pickle.loads(pickle.dumps(ox.PyFastaScanner("data/sample.fasta")))
        assert isinstance(scanner, ox.PyFastaScanner)


class TestPyFastqScanner:
    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("name", "sequence")),
        ],
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyFastqScanner("data/sample.fastq")))
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("name", "sequence", "foo"))
        scanner = pickle.loads(pickle.dumps(ox.PyFastqScanner("data/sample.fastq")))
        error = None
        try:
            scanner.scan(*input.args, **input.kwargs)
        except ValueError as e:
            error = str(e)
        finally:
            assert manifest == error

    def test_pickle(self):
        scanner = pickle.loads(pickle.dumps(ox.PyFastqScanner("data/sample.fastq")))
        assert isinstance(scanner, ox.PyFastqScanner)


class TestPyBedScanner:
    @pytest.mark.parametrize(
        "input",
        [
            *Input.permute(batch_size=[1, 2, 3, 4], bed_schema=["bed3"]),
            *Input.permute(
                batch_size=[2],
                fields=[
                    None,
                    ("chrom", "start", "end"),
                ],
                bed_schema=["bed3", "bed3+3", "bed3+6", "bed6", "bed9"],
            ),
        ],
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        key = str(input)
        scanner = pickle.loads(
            pickle.dumps(
                ox.PyBedScanner(
                    "data/sample.bed", bed_schema=input.kwargs.pop("bed_schema")
                )
            )
        )
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[key] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("nonexistent-field",))
        error = None
        try:
            scanner = pickle.loads(
                pickle.dumps(ox.PyBedScanner("data/sample.bed", bed_schema="bed9"))
            )
            schema = scanner.schema()
            stream = scanner.scan(*input.args, **input.kwargs)
            reader = pa.RecordBatchReader.from_stream(
                data=stream, schema=pa.schema(schema)
            )
            reader.read_next_batch().to_pydict()
        except ValueError as e:
            error = str(e)
        finally:
            assert manifest == error

    def test_pickle(self):
        scanner = pickle.loads(
            pickle.dumps(ox.PyBedScanner("data/sample.bed", bed_schema="bed9"))
        )
        assert isinstance(scanner, ox.PyBedScanner)


class TestPyBigBedScanner:
    @pytest.mark.parametrize(
        "input",
        [
            *Input.permute(batch_size=[1, 2, 3, 4], bed_schema=["bed3"]),
            *Input.permute(
                batch_size=[2],
                fields=[
                    None,
                    ("chrom", "start", "end"),
                    ("chrom", "start", "end", "rest"),
                ],
                bed_schema=["bed3", "bed3+3", "bed3+6", "bed6", "bed9"],
            ),
        ],
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        key = str(input)
        try:
            scanner = pickle.loads(
                pickle.dumps(
                    ox.PyBigBedScanner(
                        "data/sample.bb", schema=input.kwargs.pop("bed_schema")
                    )
                )
            )
            schema = scanner.schema()
            stream = scanner.scan(*input.args, **input.kwargs)
            reader = pa.RecordBatchReader.from_stream(
                data=stream, schema=pa.schema(schema)
            )
            batch = reader.read_next_batch()
        except Exception as e:
            assert manifest[key] == str(e)
            pass
        else:
            assert manifest[key] == batch.to_pydict()

    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("chrom", "start", "end", "chromStarts"), batch_size=2),
        ],
    )
    def test_scan_with_autosql(self, input, manifest: pytest_manifest.Manifest):
        try:
            scanner = pickle.loads(
                pickle.dumps(
                    ox.PyBigBedScanner("data/autosql-sample.bb", schema="autosql")
                )
            )
            schema = scanner.schema()
            stream = scanner.scan(*input.args, **input.kwargs)
            reader = pa.RecordBatchReader.from_stream(
                data=stream, schema=pa.schema(schema)
            )
            batch = reader.read_next_batch()
        except Exception as e:
            assert manifest[str(input)] == str(e)
            pass
        else:
            assert manifest[str(input)] == batch.to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("nonexistent-field",))
        error = None
        try:
            scanner = pickle.loads(pickle.dumps(ox.PyBigBedScanner("data/sample.bb")))
            schema = scanner.schema()
            stream = scanner.scan(*input.args, **input.kwargs)
            reader = pa.RecordBatchReader.from_stream(
                data=stream, schema=pa.schema(schema)
            )
            reader.read_next_batch().to_pydict()
        except ValueError as e:
            error = e
        finally:
            assert manifest == str(error)

    def test_pickle(self):
        scanner = pickle.loads(pickle.dumps(ox.PyBigBedScanner("data/sample.bb")))
        assert isinstance(scanner, ox.PyBigBedScanner)


class TestPyBigWigScanner:
    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("chrom", "start", "end"), batch_size=2),
        ],
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyBigWigScanner("data/sample.bw")))
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("nonexistent-field",))
        error = None
        try:
            scanner = pickle.loads(pickle.dumps(ox.PyBigWigScanner("data/sample.bw")))
            schema = scanner.schema()
            stream = scanner.scan(*input.args, **input.kwargs)
            reader = pa.RecordBatchReader.from_stream(
                data=stream, schema=pa.schema(schema)
            )
            reader.read_next_batch().to_pydict()
        except ValueError as e:
            error = e
        finally:
            assert manifest == str(error)

    def test_pickle(self):
        scanner = pickle.loads(pickle.dumps(ox.PyBigWigScanner("data/sample.bw")))
        assert isinstance(scanner, ox.PyBigWigScanner)


class TestPyGffScanner:
    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("seqid", "start", "end")),
        ],
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyGffScanner("data/sample.gff")))
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("seqid", "start", "end")),
        ],
    )
    def test_scan_with_attributes(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyGffScanner("data/sample.gff")))
        attr_defs = scanner.attribute_defs(1024)
        schema = scanner.schema(attribute_defs=attr_defs)
        stream = scanner.scan(*input.args, attribute_defs=attr_defs, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("seqid", "start", "end")),
        ],
    )
    def test_scan_sorted(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyGffScanner("data/sample.sorted.gff")))
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("seqid", "start", "end")),
        ],
    )
    def test_scan_sorted_compressed(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyGffScanner("data/sample.sorted.gff.gz", compressed=True))
        )
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    @pytest.mark.parametrize(
        "input",
        [
            *Input.permute(
                batch_size=[1, 2, 3, 4],
                fields=[None, ("seqid", "start", "end")],
                region=["chr1", "chr2"],
            ),
            Input(region="missing"),
        ],
    )
    def test_scan_query_sorted_compressed(
        self, input, manifest: pytest_manifest.Manifest
    ):
        try:
            scanner = pickle.loads(
                pickle.dumps(
                    ox.PyGffScanner("data/sample.sorted.gff.gz", compressed=True)
                )
            )
            schema = scanner.schema()
            stream = scanner.scan_query(
                *input.args,
                index="data/sample.sorted.gff.gz.tbi",
                **input.kwargs,
            )
            reader = pa.RecordBatchReader.from_stream(
                data=stream, schema=pa.schema(schema)
            )
            result = reader.read_next_batch().to_pydict()
        except Exception as e:
            result = str(e)
            pass
        assert manifest[str(input)] == result

    def test_pickle(self):
        scanner = pickle.loads(pickle.dumps(ox.PyGffScanner("data/sample.gff")))
        assert isinstance(scanner, ox.PyGffScanner)


class TestPyGtfScanner:
    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("seqid", "start", "end")),
        ],
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyGtfScanner("data/sample.gtf")))
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("seqid", "start", "end")),
        ],
    )
    def test_scan_sorted(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyGtfScanner("data/sample.sorted.gtf")))
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    @pytest.mark.parametrize(
        "input",
        [
            Input(batch_size=1),
            Input(batch_size=2),
            Input(batch_size=3),
            Input(batch_size=4),
            Input(fields=("seqid", "start", "end")),
        ],
    )
    def test_scan_sorted_compressed(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(
            pickle.dumps(ox.PyGtfScanner("data/sample.sorted.gtf.gz", compressed=True))
        )
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    @pytest.mark.parametrize(
        "input",
        [
            *Input.permute(
                batch_size=[1, 2, 3, 4],
                fields=[None, ("seqid", "start", "end")],
                region=["chr1", "chr12"],
            ),
            Input(region="missing"),
        ],
    )
    def test_scan_query_sorted_compressed(
        self, input, manifest: pytest_manifest.Manifest
    ):
        try:
            scanner = pickle.loads(
                pickle.dumps(
                    ox.PyGtfScanner("data/sample.sorted.gtf.gz", compressed=True)
                )
            )
            schema = scanner.schema()
            stream = scanner.scan_query(
                *input.args,
                index="data/sample.sorted.gtf.gz.tbi",
                **input.kwargs,
            )
            reader = pa.RecordBatchReader.from_stream(
                data=stream, schema=pa.schema(schema)
            )
            result = reader.read_next_batch().to_pydict()
        except Exception as e:
            result = str(e)
            pass
        assert manifest[str(input)] == result

    @pytest.mark.parametrize(
        "input",
        [
            Input(),
            Input(fields=("seqid", "start", "end")),
        ],
    )
    def test_schema(self, input, manifest: pytest_manifest.Manifest):
        scanner = pickle.loads(pickle.dumps(ox.PyGtfScanner("data/sample.gtf")))
        schema = scanner.schema(**input.kwargs)
        assert manifest[f"schema({str(input)})"] == schema.names

    def test_pickle(self):
        scanner = pickle.loads(pickle.dumps(ox.PyGtfScanner("data/sample.gtf")))
        assert isinstance(scanner, ox.PyGtfScanner)
