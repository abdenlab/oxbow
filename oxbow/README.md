# oxbow

The core Rust library for [oxbow](https://github.com/abdenlab/oxbow).

> **Warning**: **oxbow** is under active development. APIs are not yet stable and are subject to change.

## Installation

To use oxbow in your Rust project, add oxbow to your `Cargo.toml` or run:

```sh
cargo add oxbow
```

## Development

Ensure you have Rust installed on your system. You can install Rust using [`rustup`](https://rustup.rs/).

### Building the project

The oxbow Rust crate alone can be built using `cargo`.

```bash
cd oxbow
cargo build  # --release (for non-debug build)
```

### Linting and formatting

We use the standard Rust toolchain for linting and formatting Rust code.

[Clippy](https://doc.rust-lang.org/stable/clippy/index.html) is a Rust linter:
```bash
cargo clippy
```

The following command formats all source files of the current crate using `rustfmt`:
```bash
cargo fmt
```

### Running Tests

To run tests on Rust code, we use `cargo`:

```bash
cargo test
```
