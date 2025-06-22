---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  name: python3
---

# Marking SAM read pair duplicates with oxbow

```{code-cell} ipython3
import oxbow as ox
import polars as pl
```

## Loading the SAM File
Let's grab a sample SAM file

```{code-cell} ipython3
# Let's use oxbow to read in the SAM file as a polars dataframe
df = ox.from_sam("data/2133236").to_polars()
df.head()
```

## Compute 5' start position

```{code-cell} ipython3

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

# If the bit 0x10 is set, the read is on the reverse strand
STRAND_BIT = 0x10

def get_unclipped_5_prime_start_position(row) -> int:
    """
    Get the unclipped 5′ start position from the CIGAR string

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
        # Forward strand: 5′ end = POS - (number of leading soft-clipped bases)
        leading_soft_clips = 0
        for op, length in cigar_ops:
            if op == "S":
                leading_soft_clips += length
            else:
                break  # Stop at first non-S operation
        return pos - leading_soft_clips
    else:
        # Reverse strand: 5′ end is at POS + (aligned length) + (trailing soft-clipped bases) - 1
        aligned_length = 0
        trailing_soft_clips = 0

        # Calculate aligned length (M, =, X, D, N operations)
        for op, length in cigar_ops:
            if op in ["M", "=", "X", "D", "N"]:
                aligned_length += length

        # Find trailing soft clips (S operations at the end)
        for op, length in reversed(cigar_ops):
            if op == "S":
                trailing_soft_clips += length
            else:
                break  # Stop at first non-S operation from the end

        return pos + aligned_length + trailing_soft_clips - 1
    
df = df.with_columns(
    pl.struct(["pos", "cigar", "flag"])
    .map_elements(get_unclipped_5_prime_start_position)
    .alias("unclipped_5p_start_pos")
)

df.head()
```

## Group the reads into pairs
We assume that the qname corresponds to read pairs. We'll also be sure to include the information needed to deduplicate read pairs:
- Reference genome name
- 5' start positions
- Strand flags
- Quality scores

```{code-cell} ipython3
pairs_df = df.group_by("qname").agg(
    [
        pl.col("rname").alias("rnames"),
        pl.col("unclipped_5p_start_pos").alias("5p_positions"),
        ((pl.col("flag") & STRAND_BIT) != 0).alias("strands"),
        pl.col("qual").alias("quals"),
    ]
)
pairs_df.head()
```

## Process pair reads

Build a deduplication key and sum the quality scores for each pair (for resolution in the next step)

```{code-cell} ipython3
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


pairs_df = pairs_df.with_columns(
    pl.struct(["rnames", "5p_positions", "strands"])
    .map_elements(
        lambda s: build_dedup_key(s["rnames"], s["5p_positions"], s["strands"]),
        return_dtype=pl.String
    )
    .alias("dedup_key"),
    pl.col("quals")
    .map_elements(
        lambda qlist: sum(get_quality_score_sum(q) for q in qlist), 
        return_dtype=pl.Int64
    )
    .alias("total_quality"),
).filter(pl.col("dedup_key").is_not_null())

pairs_df[('dedup_key', 'total_quality')].head()
```

## Count the duplicates

We choose the best read pair across duplicates by the best quality score total

```{code-cell} ipython3
# Resolve duplicate pairs by sorting by total quality and taking the best pair
best_pairs_df = pairs_df.sort("total_quality", descending=True).unique(
    subset=["dedup_key"]
)

# Get the total number of duplicates
total_pair_dups = pairs_df.height - best_pairs_df.height

print(
    best_pairs_df.head(),
    "\nTotal pair duplicates:",
    total_pair_dups,
)
```
