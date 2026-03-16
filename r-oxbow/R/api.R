#' Return Arrow IPC format from a FASTQ file.
#' @export
read_fastq <- function(path, fields = "*") {
  read_fastq_impl(path, fields)
}

#' Return Arrow IPC format from a FASTA file.
#' @export
read_fasta <- function(path, regions = NULL, index = NULL, gzi = NULL, fields = "*") {
  read_fasta_impl(path, regions, index, gzi, fields)
}

#' Return Arrow IPC format from a SAM file.
#' @export
read_sam <- function(path, region = NULL, index = NULL, fields = "*", scan_rows = NULL) {
  read_sam_impl(path, region, index, fields, scan_rows)
}

#' Return Arrow IPC format from a BAM file.
#' @export
read_bam <- function(path, region = NULL, index = NULL, fields = "*", scan_rows = NULL) {
  read_bam_impl(path, region, index, fields, scan_rows)
}

#' Return Arrow IPC format from a CRAM file.
#' @export
read_cram <- function(path, reference = NULL, reference_index = NULL, region = NULL, index = NULL, fields = "*", scan_rows = NULL) {
  read_cram_impl(path, reference, reference_index, region, index, fields, scan_rows)
}

#' Return Arrow IPC format from a VCF file.
#' @export
read_vcf <- function(path, region = NULL, index = NULL, fields = "*", info_fields = "*", genotype_fields = "*", genotype_by = "sample", samples = "*", samples_nested = FALSE) {
  read_vcf_impl(path, region, index, fields, info_fields, genotype_fields, genotype_by, samples, samples_nested)
}

#' Return Arrow IPC format from a BCF file.
#' @export
read_bcf <- function(path, region = NULL, index = NULL, fields = "*", info_fields = "*", genotype_fields = "*", genotype_by = "sample", samples = "*", samples_nested = FALSE) {
  read_bcf_impl(path, region, index, fields, info_fields, genotype_fields, genotype_by, samples, samples_nested)
}

#' Return Arrow IPC format from a GTF file.
#' @export
read_gtf <- function(path, region = NULL, index = NULL, fields = "*", scan_rows = NULL) {
  read_gtf_impl(path, region, index, fields, scan_rows)
}

#' Return Arrow IPC format from a GFF file.
#' @export
read_gff <- function(path, region = NULL, index = NULL, fields = "*", scan_rows = NULL) {
  read_gff_impl(path, region, index, fields, scan_rows)
}

#' Return Arrow IPC format from a BED file.
#' @export
read_bed <- function(path, bed_schema = "bed3+", region = NULL, index = NULL, fields = "*") {
  read_bed_impl(path, bed_schema, region, index, fields)
}

#' Return Arrow IPC format from a BigWig file.
#' @export
read_bigwig <- function(path, region = NULL, fields = "*") {
  read_bigwig_impl(path, region, fields)
}

#' Return Arrow IPC format from a BigBed file.
#' @export
read_bigbed <- function(path, bed_schema = "bed3+", region = NULL, fields = "*") {
  read_bigbed_impl(path, bed_schema, region, fields)
}
