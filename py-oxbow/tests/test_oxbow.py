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
        scanner = ox.PySamScanner("data/sample.sam")
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("qname", "rname", "foo"))
        scanner = ox.PySamScanner("data/sample.sam")
        error = None
        try:
            scanner.scan(*input.args, **input.kwargs)
        except ValueError as e:
            error = str(e)
        finally:
            assert manifest == error


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
        scanner = ox.PyBamScanner("data/sample.bam")
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("qname", "rname", "foo"))
        scanner = ox.PyBamScanner("data/sample.bam")
        error = None
        try:
            scanner.scan(*input.args, **input.kwargs)
        except ValueError as e:
            error = str(e)
        finally:
            assert manifest == error


class TestPyBcfScanner:
    def test_chrom_names(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyBcfScanner("data/sample.bcf")
        assert manifest == scanner.chrom_names()

    def test_chrom_sizes(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyBcfScanner("data/sample.bcf")
        assert manifest == scanner.chrom_sizes()

    def test_info_field_defs(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyBcfScanner("data/sample.bcf")
        assert manifest == scanner.info_field_defs()

    def test_genotype_field_defs(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyBcfScanner("data/sample.bcf")
        assert manifest == scanner.genotype_field_defs()

    def test_sample_names(self):
        scanner = ox.PyBcfScanner("data/sample.bcf")
        assert 1233 == len(scanner.sample_names())

    @pytest.mark.parametrize(
        "input",
        Input.permute(batch_size=[2], samples=[["HG00096", "HG00101", "HG00103"]]),
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        scanner = ox.PyBcfScanner("data/sample.bcf")
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("name", "sequence", "foo"))
        scanner = ox.PyBcfScanner("data/sample.bcf")
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
        scanner = ox.PyBcfScanner("data/sample.bcf")
        schema = scanner.schema()
        stream = scanner.scan_query(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()


class TestPyVcfScanner:
    def test_chrom_names(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf", compressed=False)
        assert manifest == scanner.chrom_names()

    def test_chrom_names_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf.gz", compressed=True)
        assert manifest == scanner.chrom_names()

    def test_chrom_sizes(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf", compressed=False)
        assert manifest == scanner.chrom_sizes()

    def test_chrom_sizes_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf.gz", compressed=True)
        assert manifest == scanner.chrom_sizes()

    def test_info_field_names(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf", compressed=False)
        assert manifest == scanner.info_field_names()

    def test_info_field_names_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf.gz", compressed=True)
        assert manifest == scanner.info_field_names()

    def test_info_field_defs(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf", compressed=False)
        assert manifest == scanner.info_field_defs()

    def test_info_field_defs_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf.gz", compressed=True)
        assert manifest == scanner.info_field_defs()

    def test_genotype_field_names(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf", compressed=False)
        assert manifest == scanner.genotype_field_names()

    def test_genotype_field_names_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf.gz", compressed=True)
        assert manifest == scanner.genotype_field_names()

    def test_genotype_field_defs(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf", compressed=False)
        assert manifest == scanner.genotype_field_defs()

    def test_genotype_field_defs_compressed(self, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf.gz", compressed=True)
        assert manifest == scanner.genotype_field_defs()

    def test_sample_names(self):
        scanner = ox.PyVcfScanner("data/sample.vcf", compressed=False)
        assert 1233 == len(scanner.sample_names())

    def test_sample_names_compressed(self):
        scanner = ox.PyVcfScanner("data/sample.vcf.gz", compressed=True)
        assert 1233 == len(scanner.sample_names())

    @pytest.mark.parametrize(
        "input",
        Input.permute(batch_size=[2], samples=[["HG00096", "HG00101", "HG00103"]]),
    )
    def test_scan(self, input, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf", compressed=False)
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    @pytest.mark.parametrize(
        "input",
        Input.permute(batch_size=[2], samples=[["HG00096", "HG00101", "HG00103"]]),
    )
    def test_scan_compressed(self, input, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf.gz", compressed=True)
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("name", "sequence", "foo"))
        scanner = ox.PyVcfScanner("data/sample.vcf")
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
            index=["data/sample.vcf.gz.csi", "data/sample.vcf.gz.tbi"],
            batch_size=[2],
            samples=[["HG00096", "HG00101", "HG00103"]],
        ),
    )
    def test_scan_query(self, input, manifest: pytest_manifest.Manifest):
        scanner = ox.PyVcfScanner("data/sample.vcf.gz", compressed=True)
        schema = scanner.schema()
        stream = scanner.scan_query(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()


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
        scanner = ox.PyFastaScanner("data/sample.fasta")
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("name", "sequence", "foo"))
        scanner = ox.PyFastaScanner("data/sample.fasta")
        error = None
        try:
            scanner.scan(*input.args, **input.kwargs)
        except ValueError as e:
            error = str(e)
        finally:
            assert manifest == error


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
        scanner = ox.PyFastqScanner("data/sample.fastq")
        schema = scanner.schema()
        stream = scanner.scan(*input.args, **input.kwargs)
        reader = pa.RecordBatchReader.from_stream(data=stream, schema=pa.schema(schema))
        assert manifest[str(input)] == reader.read_next_batch().to_pydict()

    def test_scan_invalid_field(self, manifest):
        input = Input(fields=("name", "sequence", "foo"))
        scanner = ox.PyFastqScanner("data/sample.fastq")
        error = None
        try:
            scanner.scan(*input.args, **input.kwargs)
        except ValueError as e:
            error = str(e)
        finally:
            assert manifest == error
