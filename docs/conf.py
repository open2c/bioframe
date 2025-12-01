# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
# import sys
from datetime import datetime
from importlib.metadata import metadata

# -- Project information -----------------------------------------------------
# NOTE: If you installed your project in editable mode, this might be stale.
#       If this is the case, reinstall it to refresh the metadata
info = metadata("bioframe")
project = info["Name"]
author = "Open2C"
copyright = f"{datetime.now():%Y}, {author}."
version = info["Version"]
urls = dict(pu.split(", ") for pu in info.get_all("Project-URL"))

# The full version, including alpha/beta/rc tags
release = info["Version"]

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "autodocsumm",
    # "myst_parser",
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
]

source_suffix = [".rst", ".md"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "Thumbs.db",
    ".DS_Store",
    ".env",
    ".venv",
]

# nbsphinx_custom_formats = {
#     '.md': ['jupytext.reads', {'fmt': 'MyST'}],
# }

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

master_doc = "index"

autosummary_generate = True

# Don't include fully qualified name prefixes in autodoc
add_module_names = False

# Cache MyST (.md or .ipynb) notebook outputs if unmodified
nb_execution_mode = "cache"
nb_execution_excludepatterns = ["guide-performance.ipynb"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}
