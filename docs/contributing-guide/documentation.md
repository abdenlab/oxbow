# Documentation

Documentation requires `py-oxbow`'s `uv` environment and uses the isolated `docs` dependency group.
It is managed with the [Sphinx](https://www.sphinx-doc.org/) documentation system and uses [MyST Markdown](https://mystmd.org/) and the [Sphinx Book theme](https://sphinx-book-theme.readthedocs.io/).

To build the docs once:
```sh
uv run --group=docs sphinx-build docs docs/_build/html
```

To run a live-reload server:

```sh
uv run --group=docs sphinx-autobuild docs docs/_build/html
```