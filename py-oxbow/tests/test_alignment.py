import pytest
from pytest_manifest import Manifest
from utils import Input

import oxbow.core as ox


class TestSamFile:
    @pytest.mark.parametrize(
        "filepath",
        ["data/sample.sam", "data/malformed.sam", "data/does-not-exist.sam"],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.SamFile) as stack:
            try:
                ox.SamFile(filepath)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.SamFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [("foo",), ("foo", "bar"), ("foo", "bar", "baz"), ("*",), None],
    )
    def test_fragments(self, regions):
        for filepath in ("data/sample.bw",):
            fragments = ox.BigWigFile(filepath, regions=regions).fragments()
            assert len(fragments) == (len(regions) if regions else 1)

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            ["qname", "rname", "mapq"],
            ["qname", "rname", "foo"],
        ],
    )
    def test_batches(self, fields, manifest: Manifest):
        batches = ox.SamFile("data/sample.sam", fields=fields).batches()
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except OSError as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    def test_input_encodings(self):
        file = ox.SamFile("data/sample.sam", compressed=False, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(OSError):
            file = ox.SamFile("data/sample.sam", compressed=True, batch_size=3)
            next((file.batches()))

        file = ox.SamFile("data/sample.sam.gz", compressed=True, batch_size=3)
        assert len(next((file.batches()))) <= 3

        # file = ox.SamFile("data/sample.sam.gz", compressed=False, batch_size=3)
        # with pytest.raises(OSError):
        #     next((file.batches()))

        with pytest.raises(FileNotFoundError):
            file = ox.SamFile("doesnotexist.sam", compressed=False, batch_size=3)
            next((file.batches()))

    @pytest.mark.parametrize(
        "regions",
        [
            ["chr1"],
            ["chr1:17-32"],
            ["chr1:17-32", "chr1:30-37"],
        ],
    )
    def test_input_with_regions(self, regions):
        file = ox.SamFile(
            "data/sample.sam.gz",
            compressed=True,
            index="data/sample.sam.gz.tbi",
            regions=regions,
        )
        file.pl()

        file = ox.SamFile(
            "data/sample.sam.gz",
            compressed=True,
            index=None,  # inferred from name
            regions=regions,
        )
        file.pl()


class TestBamFile:
    @pytest.mark.parametrize(
        "filepath",
        ["data/sample.bam", "data/malformed.bam", "data/does-not-exist.bam"],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.BamFile) as stack:
            try:
                ox.BamFile(filepath, compressed=True)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f"{ox.BamFile.__name__}({Input(filepath)}, compressed=True)"
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            ["qname", "rname", "mapq"],
            ["qname", "rname", "foo"],
        ],
    )
    def test_batches(self, fields, manifest: Manifest):
        batches = ox.BamFile(
            "data/sample.bam", fields=fields, compressed=True
        ).batches()
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except OSError as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    @pytest.mark.parametrize(
        "regions",
        [["foo"], ["foo", "bar"], ["foo", "bar", "baz"], ["*"], None],
    )
    def test_fragments(self, regions):
        for filepath in ("data/sample.bam",):
            fragments = ox.BamFile(
                filepath, regions=regions, compressed=True
            ).fragments()
            assert len(fragments) == (len(regions) if regions else 1)

    def test_input_encodings(self):
        file = ox.BamFile("data/sample.bam", compressed=True, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.BamFile("data/sample.bam", compressed=False, batch_size=3)
            next((file.batches()))

        file = ox.BamFile("data/sample.ubam", compressed=False, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.BamFile("data/sample.ubam", compressed=True, batch_size=3)
            next((file.batches()))

        with pytest.raises(BaseException):
            file = ox.BamFile("doesnotexist.bam", compressed=False, batch_size=3)
            next((file.batches()))

    @pytest.mark.parametrize(
        "regions",
        [
            ["chr1"],
            ["chr1:17-32"],
            ["chr1:17-32", "chr1:30-37"],
        ],
    )
    def test_input_with_regions(self, regions):
        file = ox.BamFile(
            "data/sample.bam",
            compressed=True,
            index="data/sample.bam.bai",
            regions=regions,
        )
        file.pl()

        file = ox.BamFile(
            "data/sample.bam",
            compressed=True,
            index=None,  # inferred from name
            regions=regions,
        )
        file.pl()
