# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from __future__ import annotations

import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath(".."))

# Warning: do not change the path here. To use autodoc, you need to install the
# package first.

# -- Project information -----------------------------------------------------

project = "oxbow"
copyright = f"2023-{datetime.now().year}, Oxbow Developers"
author = "Oxbow Developers"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
]
autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# Include both markdown and rst files
source_suffix = [".rst", ".md"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "**.ipynb_checkpoints", "Thumbs.db", ".DS_Store", ".env"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
html_theme_options = {
    "logo": {
        "image_light": "_static/oxbow-logo-wide-black.svg",
        "image_dark": "_static/oxbow-logo-wide-white.svg",
    },
    "repository_url": "https://github.com/abdenlab/oxbow",
    "use_repository_button": True,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path: list[str] = ["_static"]
html_css_files: list[str] = ["style.css"]


# -- Extension configuration -------------------------------------------------
myst_enable_extensions = [
    "colon_fence",
    "deflist",
]
