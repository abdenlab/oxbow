from unittest import mock

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
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.FastaFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        ("filepath", "indexpath", "gzipath", "compressed"),
        [
            ("data/sample.fasta", "data/sample.fai", None, False),
            ("data/sample.fasta.gz", "data/sample.fai", "data/sample.fasta.gzi", True),
        ],
    )
    def test_init_regions_callstack(
        self, filepath, indexpath, gzipath, compressed, wiretap, manifest: Manifest
    ):
        inputs = Input.permute(
            [filepath],
            index=[indexpath, None],
            gzi=[gzipath, None],
            regions=[["seq1:10-20", "seq10", "seq2:1-30"]],
            compressed=[compressed],
        )
        for input in inputs:
            with wiretap(ox.FastaFile) as stack:
                try:
                    ox.FastaFile(*input.args, **input.kwargs)
                except BaseException:
                    pass
                finally:
                    assert (manifest[f"{ox.FastaFile.__name__}({input})"]) == "\n".join(
                        [c.serialize() for c in stack]
                    )

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

    @pytest.mark.parametrize(
        ("filepath", "indexpath", "gzipath", "compressed"),
        [
            ("data/sample.fasta", "data/sample.fai", None, False),
            ("data/sample.fasta.gz", "data/sample.fai", "data/sample.fasta.gzi", True),
        ],
    )
    def test_select_regions_callstack(
        self, filepath, indexpath, gzipath, compressed, wiretap, manifest: Manifest
    ):
        file_input = Input(filepath, compressed=compressed)
        file = ox.FastaFile(*file_input.args, **file_input.kwargs)
        select_inputs = Input.permute(
            index=[indexpath, None],
            gzi=[gzipath, None],
            regions=[["seq1:10-20", "seq10", "seq2:1-30"]],
        )
        for select_input in select_inputs:
            with wiretap(ox.FastaFile) as stack:
                try:
                    file.select(*select_input.args, **select_input.kwargs)
                except BaseException:
                    pass
                finally:
                    assert (
                        manifest[
                            f"{ox.FastaFile.__name__}({file_input}).select({select_input})"
                        ]
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
    @mock.patch("oxbow._core.base.pa.RecordBatchReader")
    def test_fragments(self, _, regions, mocker):
        expected_count = 1
        ds = ox.FastaFile("data/sample.fasta", index="data/sample.fai", regions=regions)
        mock_scanner_type = mocker.patch.object(ds, "_scanner_type")
        mock_scanner_type.return_value.scan_query.side_effect = regions
        fragments = ds.fragments()
        for fragment in fragments:
            fragment.iter_batches()
        assert len(fragments) == expected_count
        if regions is not None and regions != ("*",):
            assert (
                mock_scanner_type.return_value.scan_query.call_count == expected_count
            )
            for call in mock_scanner_type.return_value.scan_query.mock_calls:
                assert call.kwargs["regions"] == regions

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
