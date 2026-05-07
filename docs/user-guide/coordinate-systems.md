# Coordinate conventions

In practice, there are two commonly used conventions for the representation of genomic intervals by the numerical coordinates of their bounding bases and the semantics of inclusion or exclusion. These conventions are often informally referred to as coordinate "systems". Mixing conventions is one of the most common sources of off-by-one bugs in bioinformatics because genomic file formats and tools use or expect different conventions.

Oxbow lets you control the coordinate convention for the **output** in Arrow batches and lets you specify the interpretation of **input** query ranges, independently of the format on disk.

## Notation

SAM, VCF, and GFF use **1-based, fully-closed** intervals; BED, BigBed, and BigWig use **0-based, half-open** intervals. To make matters worse, some binary formats like BAM and BCF use 0-based representations internally while most tools surface those values as 1-based (but some don't). 

A helpful observation is that, numerically, the two conventions differ only in the encoding of the `start` coordinate. While conceptually muddled, a very useful _mnemonic_ terminology for distinguishing them is **"0-based start, 1-based end"** and **"1-based start, 1-based end"**.

Oxbow uses a compact two-character notation for coordinate conventions based on the mnemonic terminology, where the first character is the base of the start coordinate and the second is the "base" of the end coordinate:

| Code   | Informal name  | Format name | Bracket notation | Native to |
|--------|----------------|-------------|------------------|-----------|
| `"01"` | 0-based start, 1-based end | 0-based, half-open | `[start, end)` | BED, BigBed, BigWig |
| `"11"` | 1-based start, 1-based end | 1-based, closed    | `[start, end]` | SAM, BAM, CRAM, VCF, BCF, GFF, GTF † |

<small>† as returned by htslib-based tools and by noodles</small>

Oxbow defaults to each format's native convention, so you only need to think about coordinate conventions if you want output that differs from the format's native convention — for example, normalizing every source to `"01"` before joining alignment records against BED features.

## 1. Output via the `coords` argument

Every Python data source factory accepts a `coords` keyword argument:

```python
import oxbow as ox

# Native: BAM positions are emitted 1-based.
ds = ox.from_bam("data/sample.bam")
ds.pl().select("rname", "pos", "end").head()

# Coerce BAM positions to 0-based half-open to match BED tracks.
ds = ox.from_bam("data/sample.bam", coords="01")
ds.pl().select("rname", "pos", "end").head()
```

```python
# Native: BED positions are emitted 0-based.
ds = ox.from_bed("data/sample.bed")

# Coerce BED positions to 1-based closed to match SAM/VCF.
ds = ox.from_bed("data/sample.bed", coords="11")
```

Only the **start** column changes; end coordinates are the same in either convention.

## 2. Input via region queries

Representations of query regions also need to be interpreted according to a convention. Oxbow accepts two notations: implicit UCSC-style notation and explicit bracket notation.

### Implicit notation: use the DataSource's convention

Familiar `chr1:10000-20000` style, with optional `,` or `_` thousands separators:

```python
ds.regions("chr1:10,000-20,000")
ds.regions("chr1:10_000-20_000")
ds.regions("chr1")               # whole chromosome
```

This notation is **ambiguous** — `chr1:10000-20000` could mean either convention depending on context. An oxbow DataSource interprets an implicit region string according to **its own `coords` setting**. 

So `from_bam("...", coords="01").regions("chr1:10000-20000")` treats the region as 0-based half-open, matching the output it will produce.

### Bracket notation: explicit, self-describing

If you want the region to mean the same thing regardless of context, use bracket notation. The brackets carry the coordinate convention in the string itself:

```python
ds.regions("chr1:[10000,20000)")   # 0-based half-open
ds.regions("chr1:[10001,20000]")   # 1-based closed (same interval as above)
```

Bracket notation **overrides** any `coords` setting. Only `_` is accepted as a thousands separator in this form (since `,` separates start from end):

```python
ds.regions("chr1:[10_000,20_000)")
```

This is the recommended notation when a region is constructed somewhere far from the scanner that consumes it — for example, written into a config file, a CLI flag, or a pipeline manifest — because it is unambiguous on its own.

:::{tip}
A practical convention: stick to UCSC notation when you're typing regions interactively next to the scanner that uses them, and switch to bracket notation in any code that needs to be portable across formats.
:::
