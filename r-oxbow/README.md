# r-oxbow

An R interface for [oxbow](https://github.com/abdenlab/oxbow).

## Installation

> **Warning**: **oxbow** is new and under active development. It is not yet
> ready for production as APIs are subject to change.

Before you can install this package, you need to install a working [Rust toolchain](https://rustup.rs/).

On Windows, you'll also have to add the `i686-pc-windows-gnu` and `x86_64-pc-windows-gnu` targets:

```sh
rustup target add x86_64-pc-windows-gnu
rustup target add i686-pc-windows-gnu
```

Once Rust is working, you can install this package via:

```R
remotes::install_github("abdenlab/oxbow", subdir="r-oxbow")
```

> **Warning**: Using `remotes::install_github` or `devtools::install` is currently
> broken due to an [upstream issue](https://github.com/abdenlab/oxbow/issues/8).
> To try r-oxbow in the meantime, please clone the repository and run the following:

```r
rextendr::document()
devtools::load_all(".")

library(oxbow)
```

## Usage

```R
arrow_ipc <- oxbow::read_bam("./data.bam", region="chr1:1-10000")
df <- arrow::read_arrow_ipc(arrow_ipc)
head(df)
```

## Development

Changes must be recompiled as described [here](https://extendr.github.io/rextendr/articles/package.html#compile-and-use-the-package).
