# oxbow <a href="https://github.com/abdenlab/oxbow"><img align="right" src="./docs/_static/logo.svg" height="38"></img></a>
[![DOI](https://zenodo.org/badge/641466193.svg)](https://zenodo.org/doi/10.5281/zenodo.10573864)
[![Docs](https://readthedocs.org/projects/cooler/badge/?version=latest)](https://app.readthedocs.org/projects/oxbow/)

Oxbow makes genomic data ready for high-performance analytics.

Oxbow is a genomic data I/O library that models and translates next-generation sequencing (NGS) file formats into Apache [Arrow](https://arrow.apache.org/) representations. One direct application of this unification layer is the ability to access conventional NGS files as in-memory or larger-than-memory data frames in Python, R, and more. The project is organized as a multi-package monorepo with three main components:

1. [**rs-oxbow**](https://crates.io/crates/oxbow) (`oxbow/`) - A reusable, pure Rust library providing core parsing and streaming functionality
2. [**py-oxbow**](https://pypi.org/project/oxbow/) (`py-oxbow/`) - Python bindings built with PyO3/maturin
3. **r-oxbow** (`r-oxbow/`) - R bindings built with rextendr (minimal, under development)

Data I/O is handled entirely in Rust, with rich high-level features exposed via Python and R bindings.

Read the latest [documentation](https://oxbow.readthedocs.io/).

Learn more from our 
[blog post](https://open.substack.com/pub/lifeinbytes/p/breaking-out-of-bioinformatic-data-silos?r=jue12&utm_campaign=post&utm_medium=web).

## Contributing

Want to contribute? See the [contributing guide](https://oxbow.readthedocs.io/en/latest/contributing-guide/).
