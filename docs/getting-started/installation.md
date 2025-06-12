# Installation

For most users, Oxbow is a **Python** package, designed for compatibility with Arrow-based tooling in the Python data ecosystem. However, the core of Oxbow is written in Rust and the Rust library can be used for more specialized use cases or other front-end integrations. There are currently minimal bindings for R.

## Python


::::{tab-set}

:::{tab-item} PyPI
```bash
pip install oxbow
```
:::

:::{tab-item} From Source
To build and install the library from source, you will first need a working [Rust toolchain](https://rustup.rs/) as well as [maturin](https://www.maturin.rs/).

Clone the repository from GitHub:

```bash
git clone https://github.com/abdenlab/oxbow.git
```

The Python library is located in the `py-oxbow` directory. Navigate to it and install the library using `pip`:

```bash
cd py-oxbow
pip install .
```

or `uv`:

```bash
cd py-oxbow
uv sync
```

Alternatively, you can instruct `pip` to clone and install directly from GitHub by specifying the branch name and sub-directory using the following syntax:

```bash
pip install 'git+https://github.com/abdenlab/oxbow.git@main#egg=oxbow&subdirectory=py-oxbow'
```
:::

::::

## Rust


::::{tab-set}

:::{tab-item} Crates.io
```bash
cargo install oxbow
```
:::

:::{tab-item} From Source
To build and install the library from source, you will first need a working [Rust toolchain](https://rustup.rs/).

Clone the repository from GitHub:

```bash
git clone https://github.com/abdenlab/oxbow.git
```

The Rust library is located in the `oxbow` directory. Navigate to it and build the library using `cargo`:

```bash
cd oxbow
cargo build --release
```

To install the library system-wide, run:

```bash
cargo install --path .
```

:::

::::


## R

The R bindings to oxbow are currently not available on any package registry.

::::{tab-set}

:::{tab-item} From Source

Before you can install this package, you need have a working [Rust toolchain](https://rustup.rs/).
Once Rust is working, you can install this package from GitHub via `devtools::install_github` or `remotes::install_github`.

```R
> library(devtools)
> devtools::install_github("abdenlab/oxbow", subdir="r-oxbow")
```

Or clone the repository from GitHub and install locally:

```bash
git clone https://github.com/abdenlab/oxbow.git
cd oxbow
R
```

```R
> library(devtools)
> devtools::install_local("r-oxbow")
```
:::

::::
