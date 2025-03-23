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

This is a hybrid Python/Rust project and requires
[`uv`](https://github.com/astral-sh/uv) for development, relying on
[`maturin`](https://github.com/PyO3/maturin) as a build system. Dependencies
are organized into [PEP 735](https://peps.python.org/pep-0735/)-style
_dependency groups_ in the `pyproject.toml`.

By default, the development environment enables **all** dependency groups,
ensuring that all tools and libraries required by project commands are
available locally. In CI, we selectively enable groups to avoid unnecessary
compilation or installation of unused dependencies for specific commands.

```sh
uv sync # Create `.venv` with all dependency groups
```

### Building the project

To (re)build and install a local development version of `oxbow` into your
virtual environment:

```sh
uvx maturin develop --uv # --release (for non-debug build)
```

You can then test the build interactively:

```sh
uv run python
# >>> import oxbow as ox
```

or running one of the example notebooks:

```sh
uv run jupyter lab ./notebooks/bench.ipynb
```

### Running Tests

Tests use `pytest` and require only the dev dependency group (default, no
`--group` necessary):

```sh
uv run pytest
```

### Documentation

Documentation is managed with `sphinx` and uses the isolated `docs` dependency group:

```sh
uv run --group=docs sphinx-build docs docs/_build/html       # Build docs
uv run --group=docs sphinx-autobuild docs docs/_build/html   # Live-reload server
```
