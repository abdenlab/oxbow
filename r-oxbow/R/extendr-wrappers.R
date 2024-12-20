# Generated by extendr: Do not edit by hand

# nolint start

#
# This file was created with the following call:
#   .Call("wrap__make_oxbow_wrappers", use_symbols = TRUE, package_name = "oxbow")

#' @docType package
#' @usage NULL
#' @useDynLib oxbow, .registration = TRUE
NULL

#' Return Arrow IPC format from a FASTA file.
#' @export
read_fasta <- function(path, region) .Call(wrap__read_fasta, path, region)

#' Return Arrow IPC format from a FASTQ file.
#' @export
read_fastq <- function(path) .Call(wrap__read_fastq, path)

#' Return Arrow IPC format from a BAM file.
#' @export
read_bam <- function(path, region) .Call(wrap__read_bam, path, region)

#' Return Arrow IPC format from a BAM file.
#' @export
read_bam_vpos <- function(path, cpos_lo, upos_lo, cpos_hi, upos_hi) .Call(wrap__read_bam_vpos, path, cpos_lo, upos_lo, cpos_hi, upos_hi)

#' Return Arrow IPC format from a VCF file.
#' @export
read_vcf <- function(path, region) .Call(wrap__read_vcf, path, region)

#' Return Arrow IPC format from a VCF file.
#' @export
read_vcf_vpos <- function(path, cpos_lo, upos_lo, cpos_hi, upos_hi) .Call(wrap__read_vcf_vpos, path, cpos_lo, upos_lo, cpos_hi, upos_hi)

#' Return Arrow IPC format from a BCF file.
#' @export
read_bcf <- function(path, region) .Call(wrap__read_bcf, path, region)

#' Return Arrow IPC format from a BCF file.
#' @export
read_bcf_vpos <- function(path, cpos_lo, upos_lo, cpos_hi, upos_hi) .Call(wrap__read_bcf_vpos, path, cpos_lo, upos_lo, cpos_hi, upos_hi)


# nolint end
