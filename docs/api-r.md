# R API documentation

The following functions return blobs of data in Arrow IPC format (formerly known as Feather). To read the output into a `tibble` or Arrow table and interface with `dplyr` use the R [arrow](https://arrow.apache.org/docs/r/) package. For example:

```r
library(oxbow)
library(arrow)

ipc <- oxbow::read_bam("example.bam", region="chr1:1-10000")
df <- arrow::read_feather(ipc)
```

## `read_fastq`

Return Arrow IPC format from a FASTQ file.

```r
oxbow::read_fastq(path, fields = NULL)
```

## `read_fasta`

Return Arrow IPC format from a FASTA file.

```r
oxbow::read_fasta(path, regions = NULL, index = NULL, gzi = NULL, fields = NULL)
```

## `read_sam`

Return Arrow IPC format from a SAM file.

```r
oxbow::read_sam(path, region = NULL, index = NULL, fields = NULL, scan_rows = NULL)
```

## `read_bam`

Return Arrow IPC format from a BAM file.

```r
oxbow::read_bam(path, region = NULL, index = NULL, fields = NULL, scan_rows = NULL)
```

## `read_vcf`

Return Arrow IPC format from a VCF file.

```r
oxbow::read_vcf(path, region = NULL, index = NULL, fields = NULL, info_field = NULL, genotype_fields = NULL, genotype_by = "sample")
```

## `read_bcf`

Return Arrow IPC format from a BCF file.

```r
oxbow::read_bcf(path, region = NULL, index = NULL, fields = NULL, info_fields = NULL, genotype_fields = NULL, genotype_by = "sample")
```

## `read_gtf`

Return Arrow IPC format from a GTF file.

```r
oxbow::read_gtf(path, region = NULL, index = NULL, fields = NULL, scan_rows = NULL)
```

## `read_gff`

Return Arrow IPC format from a GFF file.

```r
oxbow::read_gff(path, region = NULL, index = NULL, fields = NULL, scan_rows = NULL)
```

## `read_bed`

Return Arrow IPC format from a BED file.

```r
oxbow::read_bed(path, bed_schema, region = NULL, index = NULL, fields = NULL)
```

## `read_bigwig`

Return Arrow IPC format from a BigWig file.

```r
oxbow::read_bigwig(path, region = NULL, fields = NULL)
```

## `read_bigbed`

Return Arrow IPC format from a BigBed file.

```r
oxbow::read_bigbed(path, bed_schema = "bed3+", region = NULL, fields = NULL)
```
