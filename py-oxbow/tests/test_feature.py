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
        "filepath",
        ["data/sample.bed", "data/malformed.bed", "data/does-not-exist.bed"],
    )
    def test_select_callstack(self, filepath, wiretap, manifest: Manifest):
        file = ox.BedFile(filepath)
        with wiretap(ox.BedFile) as stack:
            try:
                file.select()
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.BedFile.__name__}({Input(filepath)}).select()"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("chrom", "start", "end"), 3
        file = ox.BedFile("data/sample.bed", fields=fields)
        with wiretap(ox.BedFile) as stack:
            try:
                list(file.fragments(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.BedFile.__name__}("data/sample.bed", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_batches_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("chrom", "start", "end"), 5
        file = ox.BedFile("data/sample.bed", fields=fields)
        with wiretap(ox.BedFile) as stack:
            try:
                list(file.batches(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.BedFile.__name__}("data/sample.bed", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [("foo",), ("foo", "bar"), ("foo", "bar", "baz"), ("*",), None],
    )
    def test_fragments(self, regions):
        fragments = ox.BedFile("data/sample.bed", regions=regions).fragments()
        assert len(fragments) == (len(regions) if regions else 1)

    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (None, 3),
            (None, None),
            (("chrom", "start", "end"), 3),
            (("nonexistent-field",), 3),
        ],
    )
    def test_batches(self, fields, batch_size, manifest: Manifest):
        batches = ox.BedFile("data/sample.bed", fields=fields).batches(
            batch_size=batch_size
        )
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except ValueError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual


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
        "filepath",
        ["data/sample.bb", "data/malformed.bb", "data/does-not-exist.bb"],
    )
    def test_select_callstack(self, filepath, wiretap, manifest: Manifest):
        file = ox.BigBedFile(filepath)
        with wiretap(ox.BigBedFile) as stack:
            try:
                file.select()
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.BigBedFile.__name__}({Input(filepath)}).select()"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("chrom", "start", "end"), 3
        file = ox.BigBedFile("data/sample.bb", fields=fields)
        with wiretap(ox.BigBedFile) as stack:
            try:
                list(file.fragments(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.BigBedFile.__name__}("data/sample.bb", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_batches_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("chrom", "start", "end"), 5
        file = ox.BigBedFile("data/sample.bb", fields=fields)
        with wiretap(ox.BigBedFile) as stack:
            try:
                list(file.batches(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.BigBedFile.__name__}("data/sample.bb", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [("foo",), ("foo", "bar"), ("foo", "bar", "baz"), ("*",), None],
    )
    def test_fragments(self, regions):
        fragments = ox.BigBedFile("data/sample.bb", regions=regions).fragments()
        assert len(fragments) == (len(regions) if regions else 1)

    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (None, 3),
            (None, None),
            (("chrom", "start", "end"), 3),
            (("nonexistent-field",), 3),
        ],
    )
    def test_batches(self, fields, batch_size, manifest: Manifest):
        batches = ox.BigBedFile("data/sample.bb", fields=fields).batches(
            batch_size=batch_size
        )
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except ValueError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual


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
        "filepath",
        ["data/sample.bw", "data/malformed.bw", "data/does-not-exist.bw"],
    )
    def test_select_callstack(self, filepath, wiretap, manifest: Manifest):
        file = ox.BigWigFile(filepath)
        with wiretap(ox.BigWigFile) as stack:
            try:
                file.select()
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.BigWigFile.__name__}({Input(filepath)}).select()"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("chrom", "start", "end"), 3
        file = ox.BigWigFile("data/sample.bw", fields=fields)
        with wiretap(ox.BigWigFile) as stack:
            try:
                list(file.fragments(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.BigWigFile.__name__}("data/sample.bw", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_batches_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("chrom", "start", "end"), 5
        file = ox.BigWigFile("data/sample.bw", fields=fields)
        with wiretap(ox.BigWigFile) as stack:
            try:
                list(file.batches(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.BigWigFile.__name__}("data/sample.bw", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [("foo",), ("foo", "bar"), ("foo", "bar", "baz"), ("*",), None],
    )
    def test_fragments(self, regions):
        fragments = ox.BigWigFile("data/sample.bw", regions=regions).fragments()
        assert len(fragments) == (len(regions) if regions else 1)

    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (None, 3),
            (None, None),
            (("chrom", "start", "end"), 3),
            (("nonexistent-field",), 3),
        ],
    )
    def test_batches(self, fields, batch_size, manifest: Manifest):
        batches = ox.BigWigFile("data/sample.bw", fields=fields).batches(
            batch_size=batch_size
        )
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except ValueError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual


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
        "filepath",
        ["data/sample.gff", "data/malformed.gff", "data/does-not-exist.gff"],
    )
    def test_select_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.GffFile) as stack:
            try:
                file = ox.GffFile(filepath)
                file.select()
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.GffFile.__name__}({Input(filepath)}).select()"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("seqid", "start", "end"), 3
        file = ox.GffFile("data/sample.gff", fields=fields)
        with wiretap(ox.GffFile) as stack:
            try:
                list(file.fragments(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.GffFile.__name__}("data/sample.gff", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_batches_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("seqid", "start", "end"), 5
        file = ox.GffFile("data/sample.gff", fields=fields)
        with wiretap(ox.GffFile) as stack:
            try:
                list(file.batches(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.GffFile.__name__}("data/sample.gff", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [("foo",), ("foo", "bar"), ("foo", "bar", "baz"), ("*",), None],
    )
    def test_fragments(self, regions):
        fragments = ox.GffFile("data/sample.gff", regions=regions).fragments()
        assert len(fragments) == (len(regions) if regions else 1)

    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (None, 3),
            (None, None),
            (("seqid", "start", "end"), 3),
            (("nonexistent-field",), 3),
        ],
    )
    def test_batches(self, fields, batch_size, manifest: Manifest):
        batches = ox.GffFile("data/sample.gff", fields=fields).batches(
            batch_size=batch_size
        )
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except ValueError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual


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
        "filepath",
        ["data/sample.gtf", "data/malformed.gtf", "data/does-not-exist.gtf"],
    )
    def test_select_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.GtfFile) as stack:
            try:
                file = ox.GtfFile(filepath)
                file.select()
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.GtfFile.__name__}({Input(filepath)}).select()"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("seqid", "start", "end"), 3
        file = ox.GtfFile("data/sample.gtf", fields=fields)
        with wiretap(ox.GtfFile) as stack:
            try:
                list(file.fragments(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.GtfFile.__name__}("data/sample.gtf", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_batches_callstack(self, wiretap, manifest: Manifest):
        fields, batch_size = ("seqid", "start", "end"), 5
        file = ox.GtfFile("data/sample.gtf", fields=fields)
        with wiretap(ox.GtfFile) as stack:
            try:
                list(file.batches(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.GtfFile.__name__}("data/sample.gtf", fields={fields}).batches(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [("foo",), ("foo", "bar"), ("foo", "bar", "baz"), ("*",), None],
    )
    def test_fragments(self, regions):
        fragments = ox.GtfFile("data/sample.gtf", regions=regions).fragments()
        assert len(fragments) == (len(regions) if regions else 1)

    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (None, 3),
            (None, None),
            (("seqid", "start", "end"), 3),
            (("nonexistent-field",), 3),
        ],
    )
    def test_batches(self, fields, batch_size, manifest: Manifest):
        batches = ox.GtfFile("data/sample.gtf", fields=fields).batches(
            batch_size=batch_size
        )
        try:
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except ValueError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual

    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (None, 3),
            (("seqid", "start", "end"), 3),
            (("nonexistent-field",), 3),
        ],
    )
    def test_select(self, fields, batch_size, manifest: Manifest):
        try:
            batches = (
                ox.GtfFile("data/sample.gtf", fields=fields)
                .select()
                .to_batches(batch_size=batch_size)
            )
            actual = {f"batch-{i:02}": b.to_pydict() for i, b in enumerate(batches)}
        except OSError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual
