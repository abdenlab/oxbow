# Installation

## Python

**From PyPI**

To install the library using `pip`, run the following command:

```bash
pip install oxbow
```

**From source**

To install the library from source, clone the repository from GitHub:

```bash
git clone https://github.com/abdenlab/oxbow.git
```

The Python library is located in the `py-oxbow` directory. To install it, ensure you have Python 3.9 or later and `pip` installed on your system. Then, navigate to the `py-oxbow` directory and install the library:

Using `pip`:

```bash
cd py-oxbow
pip install .
```

Using `uv`:

```bash
cd py-oxbow
uv sync
```

## Rust

**From crates.io**

To install the library using the Rust package manager `cargo`, run the following command:

```bash
cargo install oxbow
```

**From source**

To install the library from source, clone the repository from GitHub:

```bash
git clone https://github.com/abdenlab/oxbow.git
```

The Rust library is located in the `oxbow` directory. To build and install it, ensure you have Rust installed on your system. You can install Rust using `rustup`:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

On Windows, you'll have to add the `i686-pc-windows-gnu` and `x86_64-pc-windows-gnu` targets:

```sh
rustup target add x86_64-pc-windows-gnu
rustup target add i686-pc-windows-gnu
```

Once Rust is installed, navigate to the `oxbow` directory and build the library:

```bash
cd oxbow
cargo build --release
```

To install the library system-wide, run:

```bash
cargo install --path .
```

## R

**From source**

The R bindings to oxbow are currently not available on any package registry. Before you can install this package, you need to install a working [Rust toolchain](https://rustup.rs/). Once Rust is working, you can install this package from GitHub via `devtools::install_github` or `remotes::install_github`:

```R
> remotes::install_github("abdenlab/oxbow", subdir="r-oxbow")
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
