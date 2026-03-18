---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  name: python3
---

# Quickstart

This is a quickstart guide to using Oxbow. Oxbow lets you access potentially larger-than-memory genomic files as tabular data structures, such as data frames.

## Create a DataSource

Use the convenience function associated with your file type. The returned `DataSource` object can be used to access the data in the file.

```{code-cell} ipython3
import oxbow as ox

ds = ox.from_bam("data/sample.bam")
```

### Into data frames

If the dataset fits comfortably in memory, you can materialize it fully as a [Pandas](https://pandas.pydata.org/) or [Polars](https://pola.rs/) data frame.

```{code-cell} ipython3
ds.pd()  # or ds.to_pandas()
```

```{code-cell} ipython3
ds.pl()  # or ds.to_polars()
```

### Into lazy data structures

If the data source is very large, you can also load it into a lazy or "out-of-core" data structure, such as a Polars lazy frame or [Dask](https://www.dask.org/) data frame.

```{code-cell} ipython3
df = ds.pl(lazy=True)
df.show_graph()
```

```{code-cell} ipython3
df.head().collect()
```

Oxbow data sources can also be loaded into a [DuckDB](https://duckdb.org/) relation.

```python
import duckdb

conn = duckdb.connect(":memory:")
ds = ox.from_gtf("data/gencode.v47.annotation.gtf")
rel = ds.to_duckdb(conn)
conn.sql(
    "SELECT seqid as chrom, type, start, rel.end, strand, attributes.gene_name " \
    "FROM rel " \
    "WHERE attributes.gene_name = 'PCSK9'" \
    "LIMIT 10"
).pl()
```

:::{note}
See the [Streams and Fragments](#streams-and-fragments) section for details on building Dask data frames.
:::

## Range queries

Data sources with indexes support querying genomic ranges. This is the case for htslib formats that are compressed with the BGZF gzip variant and indexed with an appropriate companion index file (e.g., `.bai`, `.tbi`, `.csi`). The BBI formats, BigWig and BigBed, possess an internal index and support range queries without an index file.

You can specify one or more ranges to the constructor or pass them to the `regions()` method. All records overlapping the query ranges will be returned.
 
```{code-cell} ipython3
ds = ox.from_bam("data/sample.bam", index="data/sample.bam.bai")
ds = ds.regions("chr1:900000-1100000")

ds.pl()
```

If the index file exists in the same location as the source file, it is automatically detected.

```{code-cell} ipython3
ox.from_bam("data/sample.bam").regions(["chr1", "chr3"]).pl()
```

```{code-cell} ipython3
ox.from_bigwig("data/sample.bw").regions("chr21:10900000-15000000").pl()
```

:::{note}
Oxbow handles multiple ranges as separate **fragments**. For more details, see the [Streams and Fragments](#streams-and-fragments) section below. 
:::

## Column projection

Oxbow lets you select only the columns you need and will not parse the others.

```{code-cell} ipython3
import polars as pl

ox.from_bam(
    "data/sample.bam", 
    fields=["rname", "pos", "end", "mapq"],
).regions(
    "chr1"
).pl()
```

In data systems lingo, selecting columns is also known as "projection". The lazy data structures returned by an oxbow data source are able to "push down" the projection operation to oxbow to prevent full record parsing. In the following example, only the four fields passed to the polars `LazyFrame.select` method will be parsed when the output gets computed.

```{code-cell} ipython3
df = (
    ox.from_bam("data/sample.bam")
    .regions("chr1")
    .pl(lazy=True)
    .select(
        pl.col("rname").alias("chrom"),
        pl.col("pos").alias("start"),
        "end",
        "mapq"
    )
    .collect()
)
df
```

## Nested and composite fields

Oxbow can handle the complex field structures of genomics file formats because they can all be mapped to Arrow constructs like lists, arrays, and structs.

For example, fields like SAM tags, VCF info and samples, and GTF attributes are exposed as struct columns in Arrow-native libraries like Polars, which are easy and efficient to manipulate.

### SAM/BAM tags

The htslib alignment formats, SAM and BAM, have optional fields called `tags` that are defined inline, rather than in a header or manifest. These definitions, a tuple of a tag name and type code, can be provided explicitly to the data source constructor for projection.

```{code-cell} ipython3
df = (
    ox.from_bam(
        "data/sample.bam", 
        fields=None,
        tag_defs=[('MD', 'Z'), ('NM', 'C')]
    )
    .regions("chr1")
    .pl()
    .select(
        pl.col("tags").struct.unnest()
    )
)
df
```

By calling the `with_tags()` method, oxbow will scan an initial number of rows to discover tag definitions to add to the schema (determined by `scan_rows`).

```{code-cell} ipython3
df = (
    ox.from_bam("data/sample.bam")
    .with_tags()
    .regions("chr1")
    .pl()
)
df
```

```{code-cell} ipython3
df['tags'].struct.unnest().head()
```

### GTF/GFF attributes

GTF/GFF attributes are analogous to SAM tags. For GTF, the type is always `"String"`. For GFF, attributes can be `"String"` or `"Array"`, the latter materializing as a list column.

```{code-cell} ipython3
df = (
    ox.from_gff("data/sample.gff")
    .with_attributes()
    .pl()
)
df.head()
```

```{code-cell} ipython3
df['attributes'].struct.unnest().head()
```

:::{important}
As of oxbow v0.7, alignment file tag definitions and annotation file attribute definitions are no longer auto-discovered by default---this behavior is **opt-in**. Use the {py:meth}`~oxbow.core.BamFile.with_tags` or {py:meth}`~oxbow.core.GffFile.with_attributes` methods, respectively, to discover or specify tag/attribute definitions.
:::


### VCF/BCF info fields

For the htslib variant call formats, VCF and BCF, the subfields of the `INFO` field are defined in the VCF header, so they do not need to be discovered by sniffing rows and you do not need to specify types.

By default, all info fields are parsed (`info_fields="*"`). You can project any subset or ignore them entirely by setting the `info_fields` argument to `None`.

```{code-cell} ipython3
(
    ox.from_vcf(
        "data/sample.vcf.gz",
        info_fields=None,
    )
    .pl()
).head()
```

```{code-cell} ipython3
df = (
    ox.from_vcf(
        "data/sample.vcf.gz",
        info_fields=["TYPE", "snpeff.Effect", "snpeff.Gene_Name", "snpeff.Transcript_BioType"],
    )
    .pl()
)
df.head()
```

```{code-cell} ipython3
df.unnest("info").head()
```

### VCF/BCF sample genotype data

For the htslib variant call formats, each variant call record is associated with an arbitrary number of so-called `FORMAT` fields that provide genotype-related information for each sample. Like `INFO`, these fields are defined in the header.

Using the `samples` and `genotype_fields` arguments, you can project any subset of samples as separate struct columns and project any subset of their associated genotype fields. Use `samples="*"` to select all samples or a list to select a subset.

```{code-cell} ipython3
df = ox.from_vcf(
    "data/sample.vcf.gz",
    info_fields=None,
    samples=['NA12891', 'NA12892'],
).pl()
df.head()
```

Each sample column is essentially a sub-dataframe of that sample's genotype fields.

```{code-cell} ipython3
df['NA12892'].struct.unnest().head()
```

:::{important}
As of oxbow v0.7, variant file sample columns are no longer projected by default---they are **opt-in**. We recommend using the {py:meth}`~oxbow.core.VcfFile.with_samples` API, below, to do this.
:::

The recommended approach to project sample genotype data is to use the `with_samples()` method. Declaring samples this way further nests all sample-related data in a single "samples" struct column for convenience.

```{code-cell} ipython3
df = (
    ox.from_vcf(
        "data/sample.vcf.gz",
        info_fields=None,
    )
    .with_samples()
    .pl()
)
df.head()
```

```{code-cell} ipython3
df.unnest("samples").head()
```

You can also customize _how_ sample genotype data are nested by using the `group_by` argument to `with_samples()`. By default (`group_by="sample"`), the columns are grouped first by sample name, then by genotype field name. By setting `group_by="field"`, you can swap the nesting order to group columns first by genotype field name, then by sample name.

```{code-cell} ipython3
df = (
    ox.from_vcf(
        "data/sample.vcf.gz",
        info_fields=None,
    )
    .with_samples(
        ['NA12891', 'NA12892'],
        genotype_fields=['AD', 'DP', 'GQ', 'PL', 'TP'],
        group_by="field",
    )    
).pl()
df.head()
```

```{code-cell} ipython3
df.unnest("samples").head()
```

In this case, each genotype field column is a data series containing the values of that field associated with each of the samples.

```{code-cell} ipython3
df.unnest("samples")['DP'].struct.unnest().head()
```

### BED schemas

Oxbow understands _BEDn+m_ schema specifiers to interpret the contents of BED files.

```{code-cell} ipython3
ox.from_bed("data/sample.bed", bed_schema="bed3+").pl().head()
```

```{code-cell} ipython3
ox.from_bed("data/sample.bed", bed_schema="bed3+6").pl().head()
```

```{code-cell} ipython3
ox.from_bed("data/sample.bed", bed_schema="bed9").pl().head()
```

### BigBed AutoSql

BigBed records natively store genomic coordinate fields and a flat string containing the "rest" of the data (equivalent to a `bed3+` schema).

```{code-cell} ipython3
ox.from_bigbed("data/autosql-sample.bb").pl().head()
```

If a BigBed file contains [AutoSql](https://genomewiki.ucsc.edu/index.php/AutoSql) definitions of its record fields and types, Oxbow can parse them.

```{code-cell} ipython3
ox.from_bigbed("data/autosql-sample.bb", schema="autosql").pl().head()
```

### Custom BED schemas

You can impose a custom parsing interpretation---field names and types (beyond the first three fields)---on a BED or BigBed file as long as the text values in those fields are compatible with the types you impose.

Pass in a BED schema as a tuple of `(str, dict[str, str])`, representing 3-12 standard BED fields (`"bed{n}"`) + custom extended fields encoded as a dictionary of field name to type name. Types can be declared using C-style AutoSql names (`string`, `short`, `float`, `double`, etc.) or Rust integer shorthands (`i8`, `u8`, `i32`, `f32`, `f64`, etc.). Fixed and variable-length array types can be declared using `int[]`, `int[10]` (AutoSql style) or `[i32]`, `[i32; 10]` (Rust shorthand style).

```{code-cell} ipython3
(
    ox.from_bigbed(
        "data/autosql-sample.bb", 
        schema=("bed4", {"score": "double", "strand": "string"})
    )
    .pl()
    .head()
)
```

```{code-cell} ipython3
narrowpeak = (
    "bed6",
    {"fold_change": "f64", "-log10p": "f64", "-log10q": "f64", "relSummit": "i64"}
)
(
    ox.from_bed(
        "data/ENCFF758CQW.100.bed.gz", 
        bed_schema=narrowpeak,
        compression="gzip"
    )
    .pl()
    .head()
)
```

## Zoom levels

The UCSC BBI formats store multiple "zoom" or "reduction" levels. These are tables of fixed-resolution genomic bins containing summary statistics of the signal of a BigWig track track or the interval coverage depth of a BigBed track.

```{code-cell} ipython3
ds = ox.from_bigwig("data/sample.bw")
ds.zoom_levels
```

```{code-cell} ipython3
ds.zoom(ds.zoom_levels[1]).regions("chr21").pl()
```

## Remote files and file-like objects

You can pull data directly from HTTP and cloud storage URLs. If needed, paths or URLs to index files must be given explicitly.

```python
ds = ox.from_bam(
    "https://oxbow-ngs.s3.us-east-2.amazonaws.com/example.bam",
    index="https://oxbow-ngs.s3.us-east-2.amazonaws.com/example.bam.bai"
)
```

Instead of using file paths or URLs, the `source` and `index` inputs to create a data source can alternatively be **callables** that open a binary I/O stream, i.e. any Python file-like object.

```python
ds = ox.from_bam(
    lambda : open("sample.bam", "rb"),
    index=lambda : open("sample.bam.bai", "rb"),
)
```

This gives you the power to customize your own transports -- to read remote sources, diverse file system implementations, or different file encodings -- independently of oxbow itself!

Libraries like [`fsspec`](https://filesystem-spec.readthedocs.io/) or [`smart_open`](https://pypi.org/project/smart-open/) can be used for this purpose.

```python
from fsspec.implementations.cached import CachingFileSystem
from s3fs import S3FileSystem

url = "https://oxbow-ngs.s3.us-east-2.amazonaws.com/example.bam"
httpfs = CachingFileSystem(target_protocol="https")
ds = ox.from_bam(
    lambda : httpfs.open(url, "rb"),
    index=lambda : httpfs.open(url + ".bai", "rb"),
)

s3fs = S3FileSystem(anon=True)
s3_uri = "s3://oxbow-ngs/example.bam"
ds = ox.from_bam(
    lambda : s3fs.open(s3_uri, "rb"),
    index=lambda : s3fs.open(s3_uri + ".bai", "rb"),
    tag_defs=[],
)
ds.regions("chr1:82744-85000").pl()
```

(streams-and-fragments)=
## Streams and Fragments

An oxbow data source object streams data via a sequence of Arrow [`RecordBatches`](https://arrow.apache.org/docs/python/generated/pyarrow.RecordBatch.html). This stream is exposed as an iterator and you can use it to materialize each batch manually. 

```{code-cell} ipython3
ds = ox.from_bam("data/sample.bam", batch_size=100)
batch = next(ds.batches())
batch
```

```{code-cell} ipython3
pl.from_arrow(batch)
```

Data sources can be logically grouped into **fragments**. Without random access, a data source contains only a single fragment. 

```{code-cell} ipython3
ds = ox.from_bam("data/sample.bam")
ds.fragments()
```

When you register range queries, each query gets mapped to a unique fragment. Each fragment generates an independent stream of record batches.

```{code-cell} ipython3
ds = ox.from_bam("data/sample.bam").regions(["chr1", "chr3", "chrX"])
ds.fragments()
```

### Dask data frames

Dask uses a different approach than the streaming paradigm of Polars and DuckDB: it subdivides a data set into a known number of independently accessible logical **partitions**, each of which is expected to fit in memory. When you convert an Oxbow data source into a Dask data frame, oxbow maps fragments to partitions:

```{code-cell} ipython3
df = (
    ox.from_bam("data/sample.bam")
    .regions(["chr1", "chrX", "chrY"])
    .dd()  # or to_dask()
)
df
```

```{code-cell} ipython3
df.partitions[1].compute()
```
