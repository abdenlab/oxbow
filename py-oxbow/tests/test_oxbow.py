from pathlib import Path

import oxbow as ox  # Remember to build via maturin in the current env
import polars as pl

# See `../../fixtures/README.md` to download files that aren't checked into the repo
test_path = Path(__file__).resolve()
project_root = test_path.parents[2]
FIXTURES_PATH = Path(project_root / "fixtures")


class TestBam:
    bam_path = str(FIXTURES_PATH / "example.bam")

    def test_read_df(self):
        ipc = ox.read_bam(self.bam_path)
        df = pl.read_ipc(ipc)

        assert not df.is_empty()

        # Check number of columns
        assert len(df.columns) == 13

    def test_read_all(self):
        ipc = ox.read_bam(self.bam_path)
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 160_178

    def test_read_region(self):
        ipc = ox.read_bam(self.bam_path, "chr1")
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 160_178

    def test_read_region_partial(self):
        ipc = ox.read_bam(self.bam_path, "chr1:1-100000")
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 4771


class TestVcf:
    vcf_path = str(FIXTURES_PATH / "ALL.chrY.phase3_integrated_v1a.20130502.genotypes.vcf.gz")

    def test_read_df(self):
        ipc = ox.read_vcf(self.vcf_path)
        df = pl.read_ipc(ipc)

        assert not df.is_empty()

        # Check number of columns
        assert len(df.columns) == 9

    def test_read_all(self):
        ipc = ox.read_vcf(self.vcf_path)
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 62_042

    def test_read_region(self):
        ipc = ox.read_vcf(self.vcf_path, "Y")
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 62_042

    def test_read_region_partial(self):
        ipc = ox.read_vcf(self.vcf_path, "Y:8028497-17629059")
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 27_947


class TestBcf:
    bcf_path = str(FIXTURES_PATH / "ALL.chrY.phase3_shapeit2_mvncall_integrated.20130502.genotypes.bcf")

    def test_read_df(self):
        ipc = ox.read_bcf(self.bcf_path)
        df = pl.read_ipc(ipc)

        assert not df.is_empty()

        # Check number of columns
        assert len(df.columns) == 9

    def test_read_all(self):
        ipc = ox.read_bcf(self.bcf_path)
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 62_042

    def test_read_region(self):
        ipc = ox.read_bcf(self.bcf_path, "Y")
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 62_042

    def test_read_region_partial(self):
        ipc = ox.read_bcf(self.bcf_path, "Y:8028497-17629059")
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 27_947


class TestBigWig:
    bigwig_path = str(FIXTURES_PATH / "valid.bigWig")

    def test_read_df(self):
        ipc = ox.read_bigwig(self.bigwig_path)
        df = pl.read_ipc(ipc)

        assert not df.is_empty()

        # Check number of columns
        assert len(df.columns) == 4

    def test_read_all(self):
        ipc = ox.read_bigwig(self.bigwig_path)
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 100_000

    def test_read_region(self):
        ipc = ox.read_bigwig(self.bigwig_path, "chr17")
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 100_000

    def test_read_region_partial(self):
        ipc = ox.read_bigwig(self.bigwig_path, "chr17:59000-60000")
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 4


class TestBigBed:
    bigbed_path = str(FIXTURES_PATH / "small.bigBed")

    def test_read_df(self):
        ipc = ox.read_bigbed(self.bigbed_path)
        df = pl.read_ipc(ipc)

        assert not df.is_empty()

        # Count number of columns
        assert len(df.columns) == 11

    def test_read_all(self):
        ipc = ox.read_bigbed(self.bigbed_path)
        df = pl.read_ipc(ipc)

        # Count number of rows
        assert len(df) == 27


class TestGff:
    gff_path = str(FIXTURES_PATH / "example.gff")

    def test_read_df(self):
        ipc = ox.read_gff(self.gff_path)
        df = pl.read_ipc(ipc)

        assert not df.is_empty()

        # Check number of columns
        assert len(df.columns) == 9

    def test_read_all(self):
        ipc = ox.read_gff(self.gff_path)
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 6


class TestGtf:
    gtf_path = str(FIXTURES_PATH / "example.gtf")

    def test_read_df(self):
        ipc = ox.read_gtf(self.gtf_path)
        df = pl.read_ipc(ipc)

        assert not df.is_empty()

        # Check number of columns
        assert len(df.columns) == 9

    def test_read_all(self):
        ipc = ox.read_gtf(self.gtf_path)
        df = pl.read_ipc(ipc)

        # Check number of rows
        assert len(df) == 2
