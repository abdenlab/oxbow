import pytest
from pytest_manifest import Manifest
from utils import Input
from unittest import mock

import oxbow.core as ox


class TestBamFile:
    @pytest.mark.parametrize(
        "filepath",
        ["data/sample.bam", "data/malformed.bam", "data/does-not-exist.bam"],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.BamFile) as stack:
            try:
                ox.BamFile(filepath)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.BamFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_select_callstack(self, wiretap, manifest: Manifest):
        bam_file = ox.BamFile("data/sample.bam")
        with wiretap(ox.BamFile) as stack:
            try:
                bam_file.select()
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.BamFile.__name__}(data/sample.bam).select()"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("qname", "rname", "mapq"), 3
        bam_file = ox.BamFile("data/sample.bam", fields=fields)
        with wiretap(ox.BamFile) as stack:
            try:
                list(bam_file.fragments(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.BamFile.__name__}("data/sample.bam", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_batches_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("qname", "rname", "mapq"), 3
        bam_file = ox.BamFile("data/sample.bam", fields=fields)
        with wiretap(ox.BamFile) as stack:
            try:
                list(bam_file.batches(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.BamFile.__name__}("data/sample.bam", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [("foo",), ("foo", "bar"), ("foo", "bar", "baz"), ("*",), None],
    )
    @mock.patch("oxbow._core.base.pa.RecordBatchReader")
    def test_fragments(self, _, regions, mocker):
        expected_count = 1 if regions is None else len(regions)
        ds = ox.BamFile("data/sample.bam", regions=regions)
        mock_scanner_type = mocker.patch.object(ds, "_scanner_type")
        mock_scanner_type.return_value.scan_query.side_effect = regions
        fragments = ds.fragments()
        for fragment in fragments:
            fragment.iter_batches()
        assert len(fragments) == expected_count
        if regions is not None and regions != ("*",):
            assert mock_scanner_type.return_value.scan_query.call_count == expected_count
            for call, region in zip(mock_scanner_type.return_value.scan_query.mock_calls, regions):
                assert call.kwargs["region"] == region

    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (None, 2),
            (None, 3),
            (None, 4),
            (None, None),
            (("qname", "rname", "mapq"), None),
            (("qname", "rname", "foo"), None),
        ],
    )
    def test_batches(self, fields, batch_size, manifest: Manifest):
        batches = ox.BamFile("data/sample.bam", fields=fields).batches(
            batch_size=batch_size
        )
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except ValueError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual


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

    def test_select_callstack(self, wiretap, manifest: Manifest):
        bam_file = ox.SamFile("data/sample.sam")
        with wiretap(ox.SamFile) as stack:
            try:
                bam_file.select()
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.SamFile.__name__}(data/sample.sam).select()"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("qname", "rname", "mapq"), 3
        bam_file = ox.SamFile("data/sample.sam", fields=fields)
        with wiretap(ox.SamFile) as stack:
            try:
                list(bam_file.fragments(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.SamFile.__name__}("data/sample.sam", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_batches_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("qname", "rname", "mapq"), 3
        bam_file = ox.SamFile("data/sample.sam", fields=fields)
        with wiretap(ox.SamFile) as stack:
            try:
                list(bam_file.batches(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.SamFile.__name__}("data/sample.sam", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [("foo",), ("foo", "bar"), ("foo", "bar", "baz"), ("*",), None],
    )
    @mock.patch("oxbow._core.base.pa.RecordBatchReader")
    def test_fragments(self, _, regions, mocker):
        expected_count = 1 if regions is None else len(regions)
        ds = ox.SamFile("data/sample.sam", regions=regions)
        mock_scanner_type = mocker.patch.object(ds, "_scanner_type")
        mock_scanner_type.return_value.scan_query.side_effect = regions
        fragments = ds.fragments()
        for fragment in fragments:
            fragment.iter_batches()
        assert len(fragments) == expected_count
        if regions is not None and regions != ("*",):
            assert mock_scanner_type.return_value.scan_query.call_count == expected_count
            for call, region in zip(mock_scanner_type.return_value.scan_query.mock_calls, regions):
                assert call.kwargs["region"] == region


    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (None, 2),
            (None, 3),
            (None, 4),
            (None, None),
            (("qname", "rname", "mapq"), None),
            (("qname", "rname", "foo"), None),
        ],
    )
    def test_batches(self, fields, batch_size, manifest: Manifest):
        batches = ox.SamFile("data/sample.sam", fields=fields).batches(
            batch_size=batch_size
        )
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except ValueError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual
