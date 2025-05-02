import pytest
from pytest_manifest import Manifest

import oxbow as ox

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
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.FastaFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "filepath",
        [
            "data/sample.fasta",
            "data/malformed.fasta",
            "data/does-not-exist.fasta",
        ],
    )
    def test_select_callstack(self, filepath, wiretap, manifest: Manifest):
        file = ox.FastaFile(filepath)
        with wiretap(ox.FastaFile) as stack:
            try:
                file.select()
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.FastaFile.__name__}({Input(filepath)}).select()"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("name", "sequence"), 3
        file = ox.FastaFile("data/sample.fasta", fields=fields)
        with wiretap(ox.FastaFile) as stack:
            try:
                list(file.fragments(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.FastaFile.__name__}("data/sample.fasta", fields={fields}).fragments(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_batches_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("name", "sequence"), 5
        file = ox.FastaFile("data/sample.fasta", fields=fields)
        with wiretap(ox.FastaFile) as stack:
            try:
                list(file.batches(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.FastaFile.__name__}("data/sample.fasta", fields={fields}).fragments(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [("foo",), ("foo", "bar"), ("foo", "bar", "baz"), ("*",), None],
    )
    def test_fragments(self, regions):
        fragments = ox.FastaFile("data/sample.fasta", regions=regions).fragments()
        assert len(fragments) == 1

    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (None, 3),
            (None, None),
            (("name", "sequence"), 3),
            (("nonexistent-field",), 3),
        ],
    )
    def test_batches(self, fields, batch_size, manifest: Manifest):
        batches = ox.FastaFile("data/sample.fasta", fields=fields).batches(
            batch_size=batch_size
        )
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except ValueError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual


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
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.FastqFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "filepath",
        [
            "data/sample.fastq",
            "data/malformed.fastq",
            "data/does-not-exist.fastq",
        ],
    )
    def test_select_callstack(self, filepath, wiretap, manifest: Manifest):
        file = ox.FastqFile(filepath)
        with wiretap(ox.FastqFile) as stack:
            try:
                file.select()
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.FastqFile.__name__}({Input(filepath)}).select()"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("name", "sequence"), 3
        file = ox.FastqFile("data/sample.fastq", fields=fields)
        with wiretap(ox.FastqFile) as stack:
            try:
                list(file.fragments(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.FastqFile.__name__}("data/sample.fastq", fields={fields}).fragments(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_batches_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("name", "sequence"), 5
        file = ox.FastqFile("data/sample.fastq", fields=fields)
        with wiretap(ox.FastqFile) as stack:
            try:
                list(file.batches(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.FastqFile.__name__}("data/sample.fastq", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments(self):
        fragments = ox.FastqFile("data/sample.fastq").fragments()
        assert len(fragments) == 1

    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (None, 3),
            (None, None),
            (("name", "sequence"), 3),
            (("nonexistent-field",), 3),
        ],
    )
    def test_batches(self, fields, batch_size, manifest: Manifest):
        batches = ox.FastqFile("data/sample.fastq", fields=fields).batches(
            batch_size=batch_size
        )
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except ValueError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual
