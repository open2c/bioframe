# Contributing

## General guidelines

If you haven't contributed to open-source before, we recommend you read [this excellent guide by GitHub on how to contribute to open source](https://opensource.guide/how-to-contribute). The guide is long, so you can gloss over things you're familiar with.

If you're not already familiar with it, we follow the [fork and pull model](https://help.github.com/articles/about-collaborative-development-models) on GitHub. Also, check out this recommended [git workflow](https://www.asmeurer.com/git-workflow/).

## Contributing Code

This project has a number of requirements for all code contributed.

- We follow the [PEP-8 style](https://www.python.org/dev/peps/pep-0008/) convention.
- We use [NumPy-style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html).
- It's ideal if user-facing API changes or new features have documentation added.
- It is best if all new functionality and/or bug fixes have unit tests added with each use-case.

## Setting up Your Development Environment

This project uses the [uv](https://docs.astral.sh/uv/) project manager and [uv_build](https://docs.astral.sh/uv/concepts/build-backend/) build backend.

To set up your dev environment, run:

```sh
uv sync
```

This command will include the requirements from the `dev` [dependency group](https://peps.python.org/pep-0735/) by default, which includes dependencies for development, testing, and documentation.

Alternatively, if you prefer to manage your virtual environments yourself, you can install the package for development using, for example:

```sh
python -m venv .venv
source .venv/bin/activate
pip install --group dev -e .
```

For all pull requests, linting and unit tests are automatically run using the [GitHub Actions](https://docs.github.com/en/actions) Continuous Integration service. However, you are still encouraged to run these checks locally before pushing code to a PR.

## Linting and Formatting

We use [ruff](https://docs.astral.sh/ruff/) for style checking.

```sh
uv run ruff check
```

Ruff can fix a lot of errors itself.

```sh
uv run ruff check --fix
```

Ruff includes a formatter that mimics [black](https://black.readthedocs.io/en/stable/). To automatically reformat your code, you can use `ruff format {source_file}`.

We use [pre-commit](https://github.com/pre-commit/pre-commit) to make sure the coding style is enforced. You first need to install pre-commit and the corresponding git commit hooks:

```sh
uv run pre-commit install
```

The last command installs the hooks listed in `.pre-commit-config.yaml` locally into your git repo. If you do this, the checks will run automatically before every commit. You can also manually make sure your code satisfies the coding style:

```sh
uv run pre-commit run --all-files
```

## Testing

It is best if all new functionality and/or bug fixes have unit tests added with each use-case.

We use [pytest](https://docs.pytest.org/en/latest) as our unit testing framework. Once you've configured your environment, you can just `cd` to the root of your repository and run:

```sh
uv run pytest
```

## Documentation

If a feature is stable and relatively finalized, it is time to add it to the documentation. If you are adding any private/public functions, it is best to add docstrings, to aid in reviewing code and also for the API reference.

We use [Numpy style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html>) and [Sphinx](http://www.sphinx-doc.org/en/stable) to document this library. Sphinx, in turn, uses [reStructuredText](http://www.sphinx-doc.org/en/stable/rest.html) as its markup language for adding code.

We use the [Sphinx Autosummary extension](http://www.sphinx-doc.org/en/stable/ext/autosummary.html) to generate API references. You may want to look at `docs/api-*.rst` files to see how they look and where to add new functions, classes or modules. We also use the [myst_nb extension](https://myst-nb.readthedocs.io/en/latest/) to render Jupyter notebooks in the documentation.

To build the documentation, run `sphinx-autobuild` using:

```sh
uv run sphinx-autobuild docs docs/_build/html
```

This will build the documentation and serve it on a local http server which listens for changes and automatically rebuilds.

Documentation from the `main` branch and tagged releases is automatically built and hosted on [readthedocs](https://readthedocs.org/).

## Acknowledgments

This document is based off of the [guidelines from the sparse project](https://github.com/pydata/sparse/blob/master/docs/contributing.rst).
