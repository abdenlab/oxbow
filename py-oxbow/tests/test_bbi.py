import cloudpickle
import fsspec
import pyarrow as pa
import pytest
from pytest_manifest import Manifest

import oxbow.core as ox
from tests.utils import Input


class TestBigBedFile:
    @pytest.mark.parametrize(
        "filepath",
        ["data/sample.bb", "data/malformed.bb", "data/does-not-exist.bb"],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.BigBedFile) as stack:
            try:
                ox.BigBedFile(filepath)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.BigBedFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [["foo"], ["foo", "bar"], ["foo", "bar", "baz"], ["*"], None],
    )
    def test_fragments(self, regions):
        for filepath in ("data/sample.bb",):
            fragments = ox.BigBedFile(filepath, regions=regions).fragments()
            assert len(fragments) == (len(regions) if regions else 1)

    def test_serialized_fragments(self):
        fragments = ox.BigBedFile(
            lambda: fsspec.open("data/sample.bb", mode="rb").open(),
            regions=["chr21"],
        ).fragments()

        fragments = cloudpickle.loads(cloudpickle.dumps(fragments))

        assert [f.count_rows() for f in fragments] == [100]

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            "*",
            ("chrom", "start", "end"),
            ("nonexistent-field",),
        ],
    )
    def test_batches(self, fields, manifest: Manifest):
        batches = ox.BigBedFile("data/sample.bb", fields=fields).batches()
        try:
            actual = {
                f"batch-{i:02}": pa.record_batch(b).to_pydict()
                for i, b in enumerate(batches)
            }
        except (OSError, ValueError) as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    def test_input_encodings(self):
        file = ox.BigBedFile("data/sample.bb", batch_size=3)
        assert next((file.batches())).num_rows <= 3

        with pytest.raises(BaseException):
            file = ox.BigBedFile("data/sample.bw", batch_size=3)
            next((file.batches()))

        with pytest.raises(BaseException):
            file = ox.BigBedFile("doesnotexist.bigBed", batch_size=3)
            next((file.batches()))

    @pytest.mark.parametrize(
        "regions",
        [
            ["chr21"],
            ["chr21:17878425-33832395", "chr21:46077848-46443662"],
        ],
    )
    def test_input_with_regions(self, regions):
        file = ox.BigBedFile("data/sample.bb", regions=regions)
        file.pl()

    def test_autosql_schema(self):
        """Autosql-derived schemas produce correct data for projected columns."""
        file = ox.BigBedFile("data/autosql-sample.bb", schema="autosql")
        batch = next(file.batches())

        # All autosql fields are present
        names = list(batch.schema.names)
        assert "chrom" in names
        assert "chromStarts" in names

        # Projection to a subset including a non-standard autosql field
        file = ox.BigBedFile(
            "data/autosql-sample.bb",
            schema="autosql",
            fields=["chrom", "start", "end", "chromStarts"],
        )
        batch = next(file.batches())
        assert list(batch.schema.names) == ["chrom", "start", "end", "chromStarts"]

        # chromStarts should have actual list data, not nulls
        col = batch.column("chromStarts")
        assert col[0].as_py() is not None, "chromStarts should not be null"
        assert isinstance(col[0].as_py(), list), "chromStarts should be a list"


class TestBigWigFile:
    @pytest.mark.parametrize(
        "filepath",
        ["data/sample.bw", "data/malformed.bw", "data/does-not-exist.bw"],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.BigWigFile) as stack:
            try:
                ox.BigWigFile(filepath)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.BigWigFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [["foo"], ["foo", "bar"], ["foo", "bar", "baz"], ["*"], None],
    )
    def test_fragments(self, regions):
        for filepath in ("data/sample.bw",):
            fragments = ox.BigWigFile(filepath, regions=regions).fragments()
            assert len(fragments) == (len(regions) if regions else 1)

    def test_serialized_fragments(self):
        fragments = ox.BigWigFile(
            lambda: fsspec.open("data/sample.bw", mode="rb").open(),
            regions=["chr21"],
        ).fragments()

        fragments = cloudpickle.loads(cloudpickle.dumps(fragments))

        assert [f.count_rows() for f in fragments] == [100]

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            "*",
            ("chrom", "start", "end"),
            ("nonexistent-field",),
        ],
    )
    def test_batches(self, fields, manifest: Manifest):
        batches = ox.BigWigFile("data/sample.bw", fields=fields).batches()
        try:
            actual = {
                f"batch-{i:02}": pa.record_batch(b).to_pydict()
                for i, b in enumerate(batches)
            }
        except (OSError, ValueError) as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    def test_input_encodings(self):
        file = ox.BigWigFile("data/sample.bw", batch_size=3)
        assert next((file.batches())).num_rows <= 3

        with pytest.raises(BaseException):
            file = ox.BigWigFile("data/sample.bb", batch_size=3)
            next((file.batches()))

        with pytest.raises(BaseException):
            file = ox.BigWigFile("doesnotexist.bigWig", batch_size=3)
            next((file.batches()))

    @pytest.mark.parametrize(
        "regions",
        [
            ["chr21"],
            ["chr21:17878425-33832395", "chr21:46077848-46443662"],
        ],
    )
    def test_input_with_regions(self, regions):
        ox.BigWigFile("data/sample.bw", regions=regions)


class TestBbiZoom:
    def test_creation(self):
        file = ox.BigWigFile("data/sample.bw")
        for res in file.zoom_levels:
            zoom = file.zoom(res, batch_size=3)
            next((zoom.batches()))

        file = ox.BigBedFile("data/sample.bb")
        for res in file.zoom_levels:
            zoom = file.zoom(res, batch_size=3)
            next((zoom.batches()))

    @pytest.mark.parametrize(
        "regions",
        [
            ["chr21"],
            ["chr21:17878425-33832395", "chr21:46077848-46443662"],
        ],
    )
    def test_creation_with_regions(self, regions):
        file = ox.BigWigFile("data/sample.bw")
        for res in file.zoom_levels:
            zoom = file.zoom(res, regions=regions, batch_size=3)
            next((zoom.batches()))

        file = ox.BigBedFile("data/sample.bb")
        for res in file.zoom_levels:
            zoom = file.zoom(res, regions=regions, batch_size=3)
            next((zoom.batches()))
