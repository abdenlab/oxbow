# py-oxbow

The Python interface for [oxbow](https://github.com/abdenlab/oxbow).

> **Warning**: **oxbow** is under active development. APIs are not yet stable and are subject to change.

## Installation

```sh
pip install oxbow
```

To build and install the bleeding edge from GitHub (requires Rust and maturin to build the package locally):

```sh
pip install 'git+https://github.com/abdenlab/oxbow.git@main#egg=oxbow&subdirectory=py-oxbow'
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

### Linting and formatting

We use [`ruff`](https://astral.sh/ruff) for linting and formatting Python code.

To validate the code style:
```sh
uv run ruff check
```

To format:
```sh
uv run ruff format
```

### Running Tests

Tests use `pytest` and require only the dev dependency group (default, no
`--group` necessary). Currently, tests must be run from within the `tests` directory to correctly generate manifests.

```sh
cd tests
uv run pytest
```

### Documentation

Documentation is managed with `sphinx` and uses the isolated `docs` dependency group:

```sh
uv run --group=docs sphinx-build docs docs/_build/html       # Build docs
uv run --group=docs sphinx-autobuild docs docs/_build/html   # Live-reload server
```
