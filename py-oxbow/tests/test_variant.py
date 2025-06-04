import pytest
from pytest_manifest import Manifest

import oxbow.core as ox
from tests.utils import Input


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
        "regions",
        [["foo"], ["foo", "bar"], ["foo", "bar", "baz"], ["*"], None],
    )
    def test_fragments(self, regions):
        for filepath in ("data/sample.vcf",):
            fragments = ox.VcfFile(filepath, regions=regions).fragments()
            assert len(fragments) == (len(regions) if regions else 1)

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            ("pos", "qual"),
            ("nonexistent-field",),
        ],
    )
    def test_batches(self, fields, manifest: Manifest):
        input = Input(
            "data/sample.vcf",
            fields=fields,
            genotype_fields=("GT",),
            info_fields=("DP",),
            samples=("HG00096",),
        )
        batches = ox.VcfFile(*input.args, **input.kwargs).batches()
        try:
            actual = len(list(batches))
        except OSError as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    def test_input_encodings(self):
        file = ox.VcfFile("data/sample.vcf", compressed=False, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.VcfFile("data/sample.vcf", compressed=True, batch_size=3)
            next((file.batches()))

        file = ox.VcfFile("data/sample.vcf.gz", compressed=True, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.VcfFile("data/sample.vcf.gz", compressed=False, batch_size=3)
            next((file.batches()))

        with pytest.raises(BaseException):
            file = ox.VcfFile("doesnotexist.vcf", compressed=False, batch_size=3)
            next((file.batches()))

    @pytest.mark.parametrize(
        "regions",
        [
            ["X"],
            ["X:51000000-51100000"],
            ["X:51000000-51100000", "X:21000000-21001000"],
        ],
    )
    def test_input_with_regions(self, regions):
        file = ox.VcfFile(
            "data/sample.vcf.gz",
            compressed=True,
            index="data/sample.vcf.gz.tbi",
            samples=["NA12878i", "NA12891", "NA12892"],
            regions=regions,
        )
        file.pl()

        file = ox.VcfFile(
            "data/sample.vcf.gz",
            compressed=True,
            index="data/sample.vcf.gz.csi",
            samples=["NA12878i", "NA12891", "NA12892"],
            regions=regions,
        )
        file.pl()

        file = ox.VcfFile(
            "data/sample.vcf.gz",
            compressed=True,
            index=None,  # inferred from name
            samples=["NA12878i", "NA12891", "NA12892"],
            regions=regions,
        )
        file.pl()


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
                ox.BcfFile(filepath, compressed=True)
            except BaseException:
                pass
            finally:
                assert (
                    manifest[
                        f"{ox.BcfFile.__name__}({Input(filepath)}, compressed=True)"
                    ]
                ) == "\n".join([c.serialize() for c in stack])

    @pytest.mark.parametrize(
        "regions",
        [["foo"], ["foo", "bar"], ["foo", "bar", "baz"], ["*"], None],
    )
    def test_fragments(self, regions):
        for filepath in ("data/sample.bcf",):
            fragments = ox.BcfFile(
                filepath, compressed=True, regions=regions
            ).fragments()
            assert len(fragments) == (len(regions) if regions else 1)

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            ("pos", "qual"),
            ("nonexistent-field",),
        ],
    )
    def test_batches(self, fields, manifest: Manifest):
        input = Input(
            "data/sample.bcf",
            compressed=True,
            fields=fields,
            genotype_fields=("GT",),
            info_fields=("DP",),
            samples=("HG00096",),
        )
        batches = ox.BcfFile(*input.args, **input.kwargs).batches()
        try:
            actual = actual = len(list(batches))
        except OSError as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    def test_input_encodings(self):
        file = ox.BcfFile("data/sample.bcf", compressed=True, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.BcfFile("data/sample.bcf", compressed=False, batch_size=3)
            next((file.batches()))

        file = ox.BcfFile("data/sample.ubcf", compressed=False, batch_size=3)
        assert len(next((file.batches()))) <= 3

        with pytest.raises(BaseException):
            file = ox.BcfFile("data/sample.ubcf", compressed=True, batch_size=3)
            next((file.batches()))

        with pytest.raises(BaseException):
            file = ox.BcfFile("doesnotexist.bcf", compressed=False, batch_size=3)
            next((file.batches()))

    @pytest.mark.parametrize(
        "regions",
        [
            ["Y"],
            ["Y:9089648-14384313"],
            ["Y:9089648-14384313", "Y:21000000-21001000"],
        ],
    )
    def test_input_with_regions(self, regions):
        file = ox.BcfFile(
            "data/sample.bcf",
            compressed=True,
            index="data/sample.bcf.csi",
            samples=["HG00096", "HG00101", "HG00103"],
            regions=regions,
        )
        file.pl()

        file = ox.BcfFile(
            "data/sample.bcf",
            compressed=True,
            index=None,  # inferred from name
            samples=["HG00096", "HG00101", "HG00103"],
            regions=regions,
        )
        file.pl()
