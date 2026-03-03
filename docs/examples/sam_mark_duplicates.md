---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  name: python3
---

# Alignment deduplication

```{code-cell} ipython3
import oxbow as ox
import polars as pl
```

## Walkthrough

To illustrate step by step, let's grab a small sample SAM file and materialize it in memory as a Polars DataFrame.

```{code-cell} ipython3
url = "https://oxbow-ngs.s3.us-east-2.amazonaws.com/Col0_C1.100k.sam"

# Let's use oxbow to read in the SAM file as a polars dataframe
df = ox.from_sam(url).to_polars()
df.head()
```

### Helper functions

```{code-cell} ipython3
# If the bit 0x10 is set, the read is on the reverse strand
STRAND_BIT = 0x10


def parse_cigar(cigar_str: str) -> list[tuple[str, int]]:
    """Parse the CIGAR string into a list of tuples (operation, length)

    Example:
        >>> parse_cigar("76M")
        [('M', 76)]
        >>> parse_cigar("10M1I65M")
        [('M', 10), ('I', 1), ('M', 65)]
    """
    result = []
    current_number = ""

    for char in cigar_str:
        if char.isdigit():
            current_number += char
        else:
            if current_number:
                result.append((char, int(current_number)))
                current_number = ""

    return result


def get_unclipped_5p_start(row) -> int:
    """
    Get the unclipped 5′ start position from the CIGAR string.

    Accounts for both soft clips (S) and hard clips (H), matching the
    reference implementation in htsjdk.

    Args:
        row: Row from the SAM file

    Returns:
        int: Unclipped 5′ start position
    """
    pos = row["pos"]
    cigar = row["cigar"]
    flag = row["flag"]

    # Parse CIGAR string
    cigar_ops = parse_cigar(cigar)

    # Check if reverse strand (bit 0x10 set)
    is_reverse = flag & STRAND_BIT

    if not is_reverse:
        # Forward strand: 5′ end = POS - (number of leading soft/hard-clipped bases)
        leading_clips = 0
        for op, length in cigar_ops:
            if op in ("S", "H"):
                leading_clips += length
            else:
                break  # Stop at first non-clip operation
        return pos - leading_clips
    else:
        # Reverse strand: 5′ end is at POS + (aligned length) + (trailing soft/hard-clipped bases) - 1
        aligned_length = 0
        trailing_clips = 0

        # Calculate aligned length (M, =, X, D, N operations)
        for op, length in cigar_ops:
            if op in ["M", "=", "X", "D", "N"]:
                aligned_length += length

        # Find trailing soft/hard clips (S or H operations at the end)
        for op, length in reversed(cigar_ops):
            if op in ("S", "H"):
                trailing_clips += length
            else:
                break  # Stop at first non-clip operation from the end

        return pos + aligned_length + trailing_clips - 1


def get_quality_score_sum(qual_str):
    """Calculate the sum of quality scores from a string of quality scores"""
    return sum(ord(c) - 33 for c in qual_str if c != " ")


def build_dedup_key(rnames, positions, strands):
    """Makes a dedup key for a read pair"""
    items = sorted(zip(rnames, positions, strands))
    if len(items) < 2:
        print(f"WARNING: read is missing pair: {items}")
        return None
    return f"{items[0][0]}:{items[0][1]}:{items[0][2]}__{items[1][0]}:{items[1][1]}:{items[1][2]}"
```

### Compute derived fields

```{code-cell} ipython3
df = df.with_columns(
    pl.struct(["pos", "cigar", "flag"])
    .map_elements(get_unclipped_5p_start, return_dtype=pl.Int64)
    .alias("5p_start"),

    pl.when((pl.col("flag") & STRAND_BIT) == 0)
    .then(pl.lit("+"))
    .otherwise(pl.lit("-"))
    .alias("strand"),

    pl.col("qual").map_elements(get_quality_score_sum, return_dtype=pl.Int64)
    .alias("total_quality")
)

df.head()
```

### Group the reads into pairs

We assume that the qname corresponds to read pairs. We group by qname and carry the original alignment records through along with the fields needed for deduplication.

```{code-cell} ipython3
pairs_df = df.group_by("qname").agg(
    [
        pl.col("rname").alias("rnames"),
        pl.col("5p_start").alias("5p_starts"),
        pl.col("strand").alias("strands"),
        pl.col("total_quality").alias("total_qualities"),
        pl.struct("*").alias("alignments"),
    ]
)
pairs_df.head()
```

### Build deduplication keys

Build a deduplication key for each read pair and filter out unpaired reads.

```{code-cell} ipython3
pairs_df = pairs_df.with_columns(
    pl.struct(["rnames", "5p_starts", "strands"])
    .map_elements(
        lambda s: build_dedup_key(s["rnames"], s["5p_starts"], s["strands"]),
        return_dtype=pl.String
    )
    .alias("dedup_key"),
).filter(pl.col("dedup_key").is_not_null())

pairs_df[("dedup_key",)].head()
```

### Resolve duplicates

We choose the best read pair across duplicates by the highest total quality score. Sorting by dedup_key first minimizes the shuffle when the data is already coordinate-sorted.

```{code-cell} ipython3
best_pairs_df = pairs_df.sort(
    ["dedup_key", "total_qualities"], descending=[False, True]
).unique(
    subset=["dedup_key"]
)

# Get the total number of duplicates
total_pair_dups = pairs_df.height - best_pairs_df.height

print("Total pair duplicates:", total_pair_dups)

best_pairs_df.head()
```

### Recover deduplicated alignments

Explode and unnest the carried alignment records to recover the original fields.

```{code-cell} ipython3
deduped_df = best_pairs_df.select(
    "qname", "alignments"
).explode("alignments").select(
    pl.col("alignments").struct.unnest()
)

deduped_df.head()
```

## Full streaming pipeline

Here's the entire deduplication pipeline chained together on a Polars LazyFrame:

```{code-cell} ipython3
ds = ox.from_sam(url)

ldf = ds.to_polars(lazy=True).with_columns(
    pl.struct(["pos", "cigar", "flag"])
    .map_elements(get_unclipped_5p_start, return_dtype=pl.Int64)
    .alias("5p_start"),

    pl.when((pl.col("flag") & STRAND_BIT) == 0)
    .then(pl.lit("+"))
    .otherwise(pl.lit("-"))
    .alias("strand"),

    pl.col("qual").map_elements(get_quality_score_sum, return_dtype=pl.Int64)
    .alias("total_quality")
).group_by("qname").agg(
    [
        pl.col("rname").alias("rnames"),
        pl.col("5p_start").alias("5p_starts"),
        pl.col("strand").alias("strands"),
        pl.col("total_quality").alias("total_qualities"),
        pl.struct(ds.schema.names).alias("alignments"),
    ]
).with_columns(
    pl.struct(["rnames", "5p_starts", "strands"])
    .map_elements(
        lambda s: build_dedup_key(s["rnames"], s["5p_starts"], s["strands"]),
        return_dtype=pl.String
    )
    .alias("dedup_key"),
).filter(
    pl.col("dedup_key").is_not_null()
).sort(
    ["dedup_key", "total_qualities"], descending=[False, True]
).unique(
    subset=["dedup_key"]
).select(
    "qname", "alignments"
).explode(
    "alignments"
).select(
    pl.col("alignments").struct.unnest()
)

ldf.show_graph()
```

Let's execute the query plan in streaming mode, writing the results to a Parquet file:

```{code-cell} ipython3
ldf.sink_parquet("data/Col0_C1.100k.dedup.pq")
```