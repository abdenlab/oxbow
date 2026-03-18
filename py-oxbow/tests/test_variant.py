import cloudpickle
import fsspec
import pyarrow as pa
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

    def test_serialized_fragments(self):
        fragments = ox.VcfFile(
            lambda: fsspec.open("data/sample.vcf.gz", mode="rb").open(),
            index=lambda: fsspec.open("data/sample.vcf.gz.tbi", mode="rb").open(),
            compressed=True,
            samples=["NA12878i", "NA12891", "NA12892"],
            regions=["X:51000000-51100000", "X:21000000-21001000"],
        ).fragments()

        fragments = cloudpickle.loads(cloudpickle.dumps(fragments))

        assert [f.count_rows() for f in fragments] == [2, 0]

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            "*",
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
            samples=("NA12878i",),
        )
        batches = ox.VcfFile(*input.args, **input.kwargs).batches()
        try:
            actual = len(list(batches))
        except Exception as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    def test_input_encodings(self):
        file = ox.VcfFile("data/sample.vcf", compressed=False, batch_size=3)
        assert next((file.batches())).num_rows <= 3

        with pytest.raises(BaseException):
            file = ox.VcfFile("data/sample.vcf", compressed=True, batch_size=3)
            next((file.batches()))

        file = ox.VcfFile("data/sample.vcf.gz", compressed=True, batch_size=3)
        assert next((file.batches())).num_rows <= 3

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

    @pytest.mark.parametrize("genotype_by", ["sample", "field"])
    def test_samples_nested_false(self, genotype_by):
        # Default: genotype columns are top-level (one column per sample or per field)
        file = ox.VcfFile(
            "data/sample.vcf",
            samples=["NA12878i", "NA12891"],
            genotype_fields=["GT"],
            genotype_by=genotype_by,
            samples_nested=False,
        )
        batch = pa.record_batch(next(file.batches()))
        assert "samples" not in batch.schema.names
        if genotype_by == "sample":
            assert "NA12878i" in batch.schema.names
            assert "NA12891" in batch.schema.names
        else:
            assert "GT" in batch.schema.names

    @pytest.mark.parametrize("genotype_by", ["sample", "field"])
    def test_samples_nested_true(self, genotype_by):
        # Nested: all genotype data wrapped under a single "samples" struct column
        file = ox.VcfFile(
            "data/sample.vcf",
            samples=["NA12878i", "NA12891"],
            genotype_fields=["GT"],
            genotype_by=genotype_by,
            samples_nested=True,
        )
        batch = pa.record_batch(next(file.batches()))
        assert "samples" in batch.schema.names
        samples_type = batch.schema.field("samples").type
        assert pa.types.is_struct(samples_type)
        if genotype_by == "sample":
            assert samples_type.get_field_index("NA12878i") >= 0
            assert samples_type.get_field_index("NA12891") >= 0
        else:
            assert samples_type.get_field_index("GT") >= 0

    def test_samples_nested_pickle_roundtrip(self):
        file = ox.VcfFile(
            "data/sample.vcf",
            samples=["NA12878i"],
            genotype_fields=["GT"],
            samples_nested=True,
        )
        file2 = cloudpickle.loads(cloudpickle.dumps(file))
        batch = pa.record_batch(next(file2.batches()))
        assert "samples" in batch.schema.names


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

    def test_serialized_fragments(self):
        fragments = ox.BcfFile(
            lambda: fsspec.open("data/sample.bcf", mode="rb").open(),
            index=lambda: fsspec.open("data/sample.bcf.csi", mode="rb").open(),
            compressed=True,
            samples=["HG00096", "HG00101", "HG00103"],
            regions=["Y:9089648-14384313", "Y:21000000-21001000"],
        ).fragments()

        fragments = cloudpickle.loads(cloudpickle.dumps(fragments))

        assert [f.count_rows() for f in fragments] == [9, 0]

    @pytest.mark.parametrize(
        "fields",
        [
            None,
            "*",
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
        except Exception as e:
            actual = str(e)

        assert manifest[f"fields={fields}"] == actual

    def test_input_encodings(self):
        file = ox.BcfFile("data/sample.bcf", compressed=True, batch_size=3)
        assert next((file.batches())).num_rows <= 3

        with pytest.raises(BaseException):
            file = ox.BcfFile("data/sample.bcf", compressed=False, batch_size=3)
            next((file.batches()))

        file = ox.BcfFile("data/sample.ubcf", compressed=False, batch_size=3)
        assert next((file.batches())).num_rows <= 3

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

    @pytest.mark.parametrize("genotype_by", ["sample", "field"])
    def test_samples_nested_false(self, genotype_by):
        file = ox.BcfFile(
            "data/sample.bcf",
            compressed=True,
            samples=["HG00096", "HG00101"],
            genotype_fields=["GT"],
            genotype_by=genotype_by,
            samples_nested=False,
        )
        batch = pa.record_batch(next(file.batches()))
        assert "samples" not in batch.schema.names
        if genotype_by == "sample":
            assert "HG00096" in batch.schema.names
            assert "HG00101" in batch.schema.names
        else:
            assert "GT" in batch.schema.names

    @pytest.mark.parametrize("genotype_by", ["sample", "field"])
    def test_samples_nested_true(self, genotype_by):
        file = ox.BcfFile(
            "data/sample.bcf",
            compressed=True,
            samples=["HG00096", "HG00101"],
            genotype_fields=["GT"],
            genotype_by=genotype_by,
            samples_nested=True,
        )
        batch = pa.record_batch(next(file.batches()))
        assert "samples" in batch.schema.names
        samples_type = batch.schema.field("samples").type
        assert pa.types.is_struct(samples_type)
        if genotype_by == "sample":
            assert samples_type.get_field_index("HG00096") >= 0
            assert samples_type.get_field_index("HG00101") >= 0
        else:
            assert samples_type.get_field_index("GT") >= 0

    def test_samples_nested_pickle_roundtrip(self):
        file = ox.BcfFile(
            "data/sample.bcf",
            compressed=True,
            samples=["HG00096"],
            genotype_fields=["GT"],
            samples_nested=True,
        )
        file2 = cloudpickle.loads(cloudpickle.dumps(file))
        batch = pa.record_batch(next(file2.batches()))
        assert "samples" in batch.schema.names
