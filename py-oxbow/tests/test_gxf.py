import pytest
from pytest_manifest import Manifest

import oxbow.core as ox
from tests.utils import Input


class TestGtfFile:
    @pytest.mark.parametrize(
        "filepath",
        ["data/sample.gtf", "data/malformed.gtf", "data/does-not-exist.gtf"],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.GtfFile) as stack:
            try:
                ox.GtfFile(filepath)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.GtfFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [["foo"], ["foo", "bar"], ["foo", "bar", "baz"], ["*"], None],
    )
    def test_fragments(self, regions):
        for filepath in ("data/sample.gtf",):
            fragments = ox.GtfFile(filepath, regions=regions).fragments()
            assert len(fragments) == (len(regions) if regions else 1)

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            ("seqid", "start", "end"),
            ("nonexistent-field",),
        ],
    )
    def test_batches(self, fields, manifest: Manifest):
        batches = ox.GtfFile("data/sample.gtf", fields=fields).batches()
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except OSError as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    def test_input_encodings(self):
        file = ox.GtfFile("data/sample.gtf", compressed=False, batch_size=3)
        assert len(next((file.batches()))) <= 3

        file = ox.GtfFile("data/sample.sorted.gtf", compressed=False, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.GtfFile("data/sample.gtf", compressed=True, batch_size=3)
            next((file.batches()))

        with pytest.raises(BaseException):
            file = ox.GtfFile("data/sample.sorted.gtf", compressed=True, batch_size=3)
            next((file.batches()))

        file = ox.GtfFile("data/sample.sorted.gtf.gz", compressed=True, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.GtfFile(
                "data/sample.sorted.gtf.gz", compressed=False, batch_size=3
            )
            next((file.batches()))

        with pytest.raises(BaseException):
            file = ox.GtfFile("doesnotexist.gtf", compressed=False, batch_size=3)
            next((file.batches()))

    @pytest.mark.parametrize(
        "regions",
        [
            ["chr1"],
            ["chr12:19000000-90000000"],
            ["chr5", "chr12:19000000-90000000"],
        ],
    )
    def test_input_with_regions(self, regions):
        file = ox.GtfFile(
            "data/sample.sorted.gtf.gz",
            compressed=True,
            index="data/sample.sorted.gtf.gz.tbi",
            regions=regions,
        )
        file.pl()

        file = ox.GtfFile(
            "data/sample.sorted.gtf.gz",
            compressed=True,
            index=None,  # inferred from name
            regions=regions,
        )
        file.pl()


class TestGffFile:
    @pytest.mark.parametrize(
        "filepath",
        ["data/sample.gff", "data/malformed.gff", "data/does-not-exist.gff"],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.GffFile) as stack:
            try:
                ox.GffFile(filepath)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.GffFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [["foo"], ["foo", "bar"], ["foo", "bar", "baz"], ["*"], None],
    )
    def test_fragments(self, regions):
        for filepath in ("data/sample.gff",):
            fragments = ox.GffFile(filepath, regions=regions).fragments()
            assert len(fragments) == (len(regions) if regions else 1)

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            ("seqid", "start", "end"),
            ("nonexistent-field",),
        ],
    )
    def test_batches(self, fields, manifest: Manifest):
        batches = ox.GffFile("data/sample.gff", fields=fields).batches()
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except OSError as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    def test_input_encodings(self):
        file = ox.GffFile("data/sample.gff", compressed=False, batch_size=3)
        assert len(next((file.batches()))) <= 3

        file = ox.GffFile("data/sample.sorted.gff", compressed=False, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.GffFile("data/sample.gff", compressed=True, batch_size=3)
            next((file.batches()))

        with pytest.raises(BaseException):
            file = ox.GffFile("data/sample.sorted.gff", compressed=True, batch_size=3)
            next((file.batches()))

        file = ox.GffFile("data/sample.sorted.gff.gz", compressed=True, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.GffFile(
                "data/sample.sorted.gff.gz", compressed=False, batch_size=3
            )
            next((file.batches()))

        with pytest.raises(BaseException):
            file = ox.GffFile("doesnotexist.gff", compressed=False, batch_size=3)
            next((file.batches()))

    @pytest.mark.parametrize(
        "regions",
        [
            ["chr2"],
            ["chr13:35000000-82000000"],
            ["chr6", "chr13:35000000-82000000"],
        ],
    )
    def test_input_with_regions(self, regions):
        file = ox.GffFile(
            "data/sample.sorted.gff.gz",
            compressed=True,
            index="data/sample.sorted.gff.gz.tbi",
            regions=regions,
        )
        file.pl()

        file = ox.GffFile(
            "data/sample.sorted.gff.gz",
            compressed=True,
            index=None,  # inferred from name
            regions=regions,
        )
        file.pl()
