# py-oxbow

A Python interface for [oxbow](https://github.com/abdenlab/oxbow).
Read specialized bioinformatic file formats into dataframes.

## Installation

> **Warning**: **oxbow** is new and under active development. It is not yet
> ready for production as APIs are subject to change.

```sh
pip install oxbow
```

## Usage

[API Documentation](https://abdenlab.org/oxbow/)

```python
import oxbow as ox

arrow_ipc = ox.read_bam("data.bam", "chr1:1-100000")

# Read into polars
import polars as pl
df = pl.read_ipc(arrow_ipc)

# Read into pandas
import io
import pyarrow.ipc
df = pyarrow.ipc.open_file(io.BytesIO(ipc)).read_pandas()
```

## Development

This project uses `maturin` and `hatch` for development, which can be installed with `pipx`.

```sh
# create a virtual env
hatch shell

# compile a development version of the package
maturin develop --release

# open the jupyter notebooks/
jupyterlab
```

A local build of `oxbow` will be added to your virtual environment.

```python
import oxbow as ox
```
