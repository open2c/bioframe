# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../"))

# autodoc_mock_imports = ["numpy", "pandas", "matplotlib", "requests"]


# -- Project information -----------------------------------------------------

project = "bioframe"
copyright = "2020, Open2C"
author = "Open2C"


# The full version, including alpha/beta/rc tags
def _read(*parts, **kwargs):
    import os

    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop("encoding", "utf-8")
    with open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text


def get_version():
    import re

    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        _read("..", "bioframe", "_version.py"),
        re.MULTILINE,
    ).group(1)
    return version


version = get_version()
# The full version, including alpha/beta/rc tags.
release = version

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # "sphinx.ext.autodoc",
    # 'sphinx.ext.doctest',
    # 'sphinx.ext.todo',
    # 'sphinx.ext.coverage',
    # 'sphinx.ext.mathjax',
    # 'sphinx.ext.ifconfig',
    "autodocsumm",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",  # 'numpydoc'
    "myst_nb",
]
# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", '**.ipynb_checkpoints']

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
jupyter_execute_notebooks = "cache"
execution_excludepatterns = ['guide-performance.ipynb']
