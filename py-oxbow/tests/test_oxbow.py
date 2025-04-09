import pyarrow as pa

import pytest
import pytest_manifest

import oxbow.oxbow as ox

from utils import Input


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
