# Format matrix

Oxbow currently supports the following NGS formats, with the encodings and access patterns described in the table below. Oxbow enables random access to most uncompressed or structured, block-compressed data sources and can use conventional format-specific indexes to support efficient genomic range queries.

| Format    | Modality | Encoding | Multiple resolutions | Compression support[^1] | Index support | Authority |
|-----------|----------|----------|----------------------|-------------------------|---------------|-----------|
| FASTQ     | Sequence | Text     | No                  | BGZF                     |               |           |
| FASTA     | Sequence | Text     | No                  | BGZF                     | FAI, FAI + GZI[^2]  |     |
| SAM       | Alignment| Text     | No                  | BGZF                     | TBI, CSI      | [3] |
| BAM       | Alignment| Binary   | No                  | BGZF                     | BAI, CSI      | [3] |
| CRAM      | Alignment| Binary   | No                  | many[^6]                 | CRAI          | [3] |
| VCF       | Variant call | Text | No                  | BGZF                     | TBI, CSI      | [3] |
| BCF       | Variant call | Binary | No                | BGZF                     | CSI           | [3] |
| GTF       | Feature  | Text     | No                  | BGZF                     | TBI, CSI      |     |
| GFF3      | Feature  | Text     | No                  | BGZF                     | TBI, CSI      | [5] |
| BED       | Feature  | Text     | No                  | BGZF                     | TBI, CSI       | [3], [4] |
| BigBed    | Feature  | Binary   | Yes                 | deflate                  | Internal index |  [4] |
| BigWig    | Feature  | Binary   | Yes                 | deflate                  | Internal index | [4] |

[^1]: BGZF is a backwards-compatible "blocked" variant of GZIP that allows for random access and indexing using one of the companion index files listed in the next column.
[^2]: FAI indexes uncompressed FASTA records. GZI is an additional file required to perform look ups in BGZF-compressed FASTA.
[^3]: GA4GH and the [htslib](https://www.htslib.org/) references and specifications.
[^4]: UCSC Genome Browser and the [Kent Source](https://github.com/ucscGenomeBrowser/kent-core).
[^5]: [sequenceontology.org](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
[^6]: CRAM is a hierarchical format, organized into "containers" and "slices", with slices being the unit of random access. Each slice consists of data "blocks" that support a number of different encodings, including lossy ones. For example, for CRAM files using reference-based compression, the exact (external) reference sequences are required to decode sequence bases from alignment records.