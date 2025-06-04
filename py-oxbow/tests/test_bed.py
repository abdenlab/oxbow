import pytest
from pytest_manifest import Manifest

import oxbow.core as ox
from tests.utils import Input


class TestBedFile:
    @pytest.mark.parametrize(
        "filepath",
        ["data/sample.bed", "data/malformed.bed", "data/does-not-exist.bed"],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.BedFile) as stack:
            try:
                ox.BedFile(filepath)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.BedFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [["foo"], ["foo", "bar"], ["foo", "bar", "baz"], ["*"], None],
    )
    def test_fragments(self, regions):
        for filepath in ("data/sample.bed",):
            fragments = ox.BedFile(filepath, regions=regions).fragments()
            assert len(fragments) == (len(regions) if regions else 1)

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            ("chrom", "start", "end"),
            ("nonexistent-field",),
        ],
    )
    def test_batches(self, fields, manifest: Manifest):
        batches = ox.BedFile("data/sample.bed", fields=fields).batches()
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except OSError as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    def test_input_encodings(self):
        file = ox.BedFile("data/sample.bed", "bed9", compressed=False, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.BedFile("data/sample.bed", "bed9", compressed=True, batch_size=3)
            next((file.batches()))

        file = ox.BedFile("data/sample.bed.gz", "bed9", compressed=True, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.BedFile(
                "data/sample.bed.gz", "bed9", compressed=False, batch_size=3
            )
            next((file.batches()))

        with pytest.raises(BaseException):
            file = ox.BedFile(
                "doesnotexist.bed", "bed9", compressed=False, batch_size=3
            )
            next((file.batches()))

    @pytest.mark.parametrize(
        "regions",
        [
            ["chr1"],
            ["chr12"],
            ["chr5", "chr12:1-100"],
        ],
    )
    def test_input_with_regions(self, regions):
        file = ox.BedFile(
            "data/sample.bed.gz",
            compressed=True,
            index="data/sample.bed.gz.tbi",
            regions=regions,
        )
        file.pl()

        file = ox.BedFile(
            "data/sample.bed.gz",
            compressed=True,
            index=None,  # inferred from name
            regions=regions,
        )
        file.pl()
