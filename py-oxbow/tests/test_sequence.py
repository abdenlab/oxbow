import pytest
from pytest_manifest import Manifest

import oxbow.core as ox
from tests.utils import Input


class TestFastaFile:
    @pytest.mark.parametrize(
        "filepath",
        [
            "data/sample.fasta",
            "data/malformed.fasta",
            "data/does-not-exist.fasta",
        ],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.FastaFile) as stack:
            try:
                ox.FastaFile(filepath)
            finally:
                assert (
                    manifest[f"{ox.FastaFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [["foo"], ["foo", "bar"], ["foo", "bar", "baz"], ["*"], None],
    )
    def test_fragments(self, regions):
        fragments = ox.FastaFile(
            "data/sample.fasta", index="data/sample.fasta.fai", regions=regions
        ).fragments()
        assert len(fragments) == 1

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            ("name", "sequence"),
            ("nonexistent-field",),
        ],
    )
    def test_batches(self, fields, manifest: Manifest):
        batches = ox.FastaFile("data/sample.fasta", fields=fields).batches()
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except OSError as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    @pytest.mark.parametrize(
        ("filepath", "indexpath", "gzipath", "compressed", "regions"),
        [
            (
                "data/sample.fasta",
                "data/sample.fasta.fai",
                None,
                False,
                ["seq1:10-20", "seq10"],
            ),
            (
                "data/sample.fasta",
                "data/sample.fasta.fai",
                None,
                False,
                ["seq1:10-20", "seq10", "seq2:1-30"],
            ),
            (
                "data/sample.fasta.gz",
                "data/sample.fasta.fai",
                "data/sample.fasta.gz.gzi",
                True,
                ["seq1:10-20", "seq10"],
            ),
            (
                "data/sample.fasta.gz",
                "data/sample.fasta.fai",
                "data/sample.fasta.gz.gzi",
                True,
                ["seq1:10-20", "seq10", "seq2:1-30"],
            ),
        ],
    )
    def test_batches_with_regions(
        self, filepath, indexpath, gzipath, compressed, regions, manifest: Manifest
    ):
        input = Input(
            filepath,
            index=indexpath,
            gzi=gzipath,
            compressed=compressed,
            regions=regions,
        )
        batches = ox.FastaFile(*input.args, **input.kwargs).batches()
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except ValueError as e:
            actual = str(e)

        assert manifest[f"{ox.FastaFile.__name__}({input}).batches()"] == actual

    def test_input_encodings(self):
        file = ox.FastaFile("data/sample.fasta", compressed=False, batch_size=3)
        assert len(next((file.batches()))) <= 3

        file = ox.FastaFile("data/sample.fasta", compressed=True, batch_size=3)
        with pytest.raises(OSError):
            next((file.batches()))

        file = ox.FastaFile("data/sample.fasta.gz", compressed=True, batch_size=3)
        assert len(next((file.batches()))) <= 3

        file = ox.FastaFile("data/sample.fasta.gz", compressed=False, batch_size=3)
        with pytest.raises(OSError):
            next((file.batches()))

        file = ox.FastaFile("doesnotexist.fasta", compressed=False, batch_size=3)
        with pytest.raises(FileNotFoundError):
            next((file.batches()))

    @pytest.mark.parametrize(
        "regions",
        [
            ["seq1"],
            ["seq1:10-20"],
            ["seq1:10-20", "seq10"],
            ["seq1:10-20", "seq10", "seq2:1-30"],
        ],
    )
    def test_input_with_regions(self, regions):
        file = ox.FastaFile(
            "data/sample.fasta",
            compressed=False,
            index="data/sample.fasta.fai",
            regions=regions,
        )
        assert len(file.pl()) == len(regions)

        file = ox.FastaFile(
            "data/sample.fasta",
            compressed=False,
            index=None,  # inferred from name
            regions=regions,
        )
        assert len(file.pl()) == len(regions)

        file = ox.FastaFile(
            "data/sample.fasta.gz",
            compressed=True,
            index="data/sample.fasta.fai",
            gzi="data/sample.fasta.gz.gzi",
            regions=regions,
        )
        assert len(file.pl()) == len(regions)

        file = ox.FastaFile(
            "data/sample.fasta.gz",
            compressed=True,
            index="data/sample.fasta.fai",
            gzi=None,  # currently not inferred from name
            regions=regions,
        )
        with pytest.raises(ValueError):
            file.pl()


class TestFastqFile:
    @pytest.mark.parametrize(
        "filepath",
        [
            "data/sample.fastq",
            "data/malformed.fastq",
            "data/does-not-exist.fastq",
        ],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.FastqFile) as stack:
            try:
                ox.FastqFile(filepath)
            finally:
                assert (
                    manifest[f"{ox.FastqFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments(self):
        fragments = ox.FastqFile("data/sample.fastq").fragments()
        assert len(fragments) == 1

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            ("name", "sequence"),
            ("nonexistent-field",),
        ],
    )
    def test_batches(self, fields, manifest: Manifest):
        batches = ox.FastqFile("data/sample.fastq", fields=fields).batches()
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except OSError as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    @pytest.mark.parametrize(
        "filepath",
        [
            "data/sample.fastq",
            "data/malformed.fastq",
            "data/does-not-exist.fastq",
        ],
    )
    def test_select(self, filepath):
        file = ox.FastqFile(filepath)
        with pytest.raises(NotImplementedError):
            file.regions("seq1:10-20")

    def test_input_encodings(self):
        file = ox.FastqFile("data/sample.fastq", compressed=False, batch_size=3)
        assert len(next((file.batches()))) <= 3

        file = ox.FastqFile("data/sample.fastq", compressed=True, batch_size=3)
        with pytest.raises(OSError):
            next((file.batches()))

        file = ox.FastqFile("data/sample.fastq.gz", compressed=True, batch_size=3)
        assert len(next((file.batches()))) <= 3

        file = ox.FastqFile("data/sample.fastq.gz", compressed=False, batch_size=3)
        with pytest.raises(OSError):
            next((file.batches()))

        file = ox.FastqFile("data/malformed.fastq", compressed=False, batch_size=3)
        with pytest.raises(OSError):
            next((file.batches()))

        file = ox.FastqFile("doesnotexist.fastq", compressed=False, batch_size=3)
        with pytest.raises(FileNotFoundError):
            next((file.batches()))
