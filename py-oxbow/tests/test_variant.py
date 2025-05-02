import pytest
from pytest_manifest import Manifest

import oxbow.core as ox

from tests.utils import Input


class TestBcfFile:
    @pytest.mark.parametrize(
        "filepath",
        [
            "data/sample.bcf",
            "data/malformed.bcf",
            "data/does-not-exist.bcf",
        ],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.BcfFile) as stack:
            try:
                ox.BcfFile(filepath)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.BcfFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "filepath",
        [
            "data/sample.bcf",
            "data/malformed.bcf",
            "data/does-not-exist.bcf",
        ],
    )
    def test_select_callstack(self, filepath, wiretap, manifest: Manifest):
        input = Input(
            filepath,
            samples=("HG00096", "HG00101", "HG00103"),
            genotype_fields=("GT",),
            info_fields=("DP",),
        )
        with wiretap(ox.BcfFile) as stack:
            try:
                file = ox.BcfFile(*input.args, **input.kwargs)
                file.select()
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.BcfFile.__name__}({Input(filepath)}).select()"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments_callstack(self, wiretap, manifest: Manifest):
        samples, batch_size = ("HG00096", "HG00101", "HG00103"), 3
        file = ox.BcfFile("data/sample.bcf", samples=samples)
        with wiretap(ox.BcfFile) as stack:
            try:
                list(file.fragments(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.BcfFile.__name__}("data/sample.bcf", samples={samples}).fragments(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_batches_callstack(self, wiretap, manifest: Manifest):
        batch_size = 3
        input = Input(
            "data/sample.bcf",
            samples=("HG00096", "HG00101", "HG00103"),
            genotype_fields=("GT",),
            info_fields=("DP",),
        )
        with wiretap(ox.BcfFile) as stack:
            try:
                file = ox.BcfFile(*input.args, **input.kwargs)
                file.batches(batch_size=batch_size)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f"{ox.BcfFile.__name__}({str(input)}).fragments(batch_size={batch_size})"
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [("foo",), ("foo", "bar"), ("foo", "bar", "baz"), ("*",), None],
    )
    def test_fragments(self, regions):
        fragments = ox.BcfFile(
            "data/sample.bcf",
            regions=regions,
            samples=("HG00096", "HG00101", "HG00103"),
        ).fragments()
        assert len(fragments) == len(regions) if regions else 1

    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (("pos", "qual"), 3),
            (("nonexistent-field",), 3),
        ],
    )
    def test_batches(self, fields, batch_size, manifest: Manifest):
        input = Input(
            "data/sample.bcf",
            fields=fields,
            genotype_fields=("GT",),
            info_fields=("DP",),
            samples=("HG00096",),
        )
        batches = ox.BcfFile(*input.args, **input.kwargs).batches(batch_size=batch_size)
        try:
            actual = actual = len(list(batches))
        except ValueError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual


class TestVcfFile:
    @pytest.mark.parametrize(
        "filepath",
        [
            "data/sample.vcf",
            "data/malformed.vcf",
            "data/does-not-exist.vcf",
        ],
    )
    def test_init_callstack(self, filepath, wiretap, manifest: Manifest):
        with wiretap(ox.VcfFile) as stack:
            try:
                ox.VcfFile(filepath)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.VcfFile.__name__}({Input(filepath)})"]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "filepath",
        [
            "data/sample.vcf",
            "data/malformed.vcf",
            "data/does-not-exist.vcf",
        ],
    )
    def test_select_callstack(self, filepath, wiretap, manifest: Manifest):
        input = Input(
            filepath,
            samples=("HG00096", "HG00101", "HG00103"),
            genotype_fields=("GT",),
            info_fields=("DP",),
        )
        with wiretap(ox.VcfFile) as stack:
            try:
                file = ox.VcfFile(*input.args, **input.kwargs)
                file.select()
            except BaseException:
                pass
            finally:
                assert (
                    manifest[f"{ox.VcfFile.__name__}({input}).select()"]
                ) == "\n".join([c.serialize() for c in stack])

    def test_fragments_callstack(self, wiretap, manifest: Manifest):
        samples, batch_size = ("HG00096", "HG00101", "HG00103"), 3
        file = ox.VcfFile("data/sample.vcf", samples=samples)
        with wiretap(ox.VcfFile) as stack:
            try:
                list(file.fragments(batch_size=batch_size))
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f'{ox.VcfFile.__name__}("data/sample.vcf", samples={samples}).fragments(batch_size={batch_size})'
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    def test_batches_callstack(self, wiretap, manifest: Manifest):
        batch_size = 3
        input = Input(
            "data/sample.vcf",
            samples=("HG00096", "HG00101", "HG00103"),
            genotype_fields=("GT",),
            info_fields=("DP",),
        )
        with wiretap(ox.VcfFile) as stack:
            try:
                file = ox.VcfFile(*input.args, **input.kwargs)
                file.batches(batch_size=batch_size)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f"{ox.VcfFile.__name__}({str(input)}).fragments(batch_size={batch_size})"
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [("foo",), ("foo", "bar"), ("foo", "bar", "baz"), ("*",), None],
    )
    def test_fragments(self, regions):
        fragments = ox.VcfFile(
            "data/sample.vcf",
            regions=regions,
            samples=("HG00096", "HG00101", "HG00103"),
        ).fragments()
        assert len(fragments) == len(regions) if regions else 1

    @pytest.mark.parametrize(
        ("fields", "batch_size"),
        [
            (None, 1),
            (("pos", "qual"), 3),
            (("nonexistent-field",), 3),
        ],
    )
    def test_batches(self, fields, batch_size, manifest: Manifest):
        input = Input(
            "data/sample.vcf",
            fields=fields,
            genotype_fields=("GT",),
            info_fields=("DP",),
            samples=("HG00096",),
        )
        batches = ox.VcfFile(*input.args, **input.kwargs).batches(batch_size=batch_size)
        try:
            actual = len(list(batches))
        except ValueError as e:
            actual = str(e)

        assert manifest[f"fields={fields}, batch_size={batch_size}"] == actual
