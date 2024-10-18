# Contributing


## General guidelines

If you haven't contributed to open-source before, we recommend you read [this excellent guide by GitHub on how to contribute to open source](https://opensource.guide/how-to-contribute). The guide is long, so you can gloss over things you're familiar with.

If you're not already familiar with it, we follow the [fork and pull model](https://help.github.com/articles/about-collaborative-development-models) on GitHub. Also, check out this recommended [git workflow](https://www.asmeurer.com/git-workflow/).


## Contributing Code

This project has a number of requirements for all code contributed.

* We follow the [PEP-8 style](https://www.python.org/dev/peps/pep-0008/) convention.
* We use [NumPy-style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html).
* It's ideal if user-facing API changes or new features have documentation added.
* It is best if all new functionality and/or bug fixes have unit tests added with each use-case.


## Setting up Your Development Environment

This project uses the [hatch](https://hatch.pypa.io/latest/) project manager and build system. We recommend you install `hatch` as a global isolated application using [pipx](https://pipx.pypa.io/stable/). See other installation options [here](https://hatch.pypa.io/latest/install/).

```sh
pipx install hatch
```

> [!NOTE]
> Many custom command shortcuts are accessible through hatch (and shown below). See `tool.hatch.envs.default.scripts` in our project's `pyproject.toml` configuration file.

After forking and cloning the repository, you can create an isolated Python development environment and install the package in "editable" (i.e. development) mode as follows:

```sh
git clone https://github.com/open2c/bioframe.git
cd bioframe
hatch shell
```

The first time you run `hatch shell` the environment will be created and activated, and the package will be installed. In future sessions, running `hatch shell` will reactivate your development environment.

> [!TIP]
> If you prefer to store your virtual environments in your working directory (like classic virtualenvs) rather than in a centralized location (similar to conda), configure hatch as follows:
>
> ```sh
> hatch config set dirs.env.virtual .venv
> ```
>
> This will make hatch set up its environments within the current working directory under `.venv`.

Alternatively, if you prefer to manage your virtual environments yourself, you can install the package for development using, for example:

```sh
python -m venv .venv
source .venv/bin/activate
pip install -e '.[dev,test,docs]'
```

For all pull requests, linting and unit tests are automatically run using the [GitHub Actions](https://docs.github.com/en/actions) Continuous Integration service. However, you are still encouraged to run these checks locally before pushing code to a PR.

## Linting and Formatting

We use [ruff](https://docs.astral.sh/ruff/) for style checking. Run `ruff check .` or:

```sh
hatch run lint
```

Ruff can fix a lot of errors itself. Run `ruff check --fix .` or:

```sh
hatch run fix
```

Ruff includes a formatter that mimics [black](https://black.readthedocs.io/en/stable/). To automatically reformat your code, you can use `ruff format {source_file}`.

We use [pre-commit](https://github.com/pre-commit/pre-commit) to make sure the coding style is enforced. You first need to install pre-commit and the corresponding git commit hooks:

```sh
pip install pre-commit
pre-commit install
```

The last command installs the hooks listed in `.pre-commit-config.yaml` locally into your git repo. If you do this, the checks will run automatically before every commit. You can also manually make sure your code satisfies the coding style:

```sh
pre-commit run --all-files
```

## Testing

It is best if all new functionality and/or bug fixes have unit tests added with each use-case.

We use [pytest](https://docs.pytest.org/en/latest) as our unit testing framework. Once you've configured your environment, you can just `cd` to the root of your repository and run `pytest` or:

```sh
hatch run test
```

## Documentation

If a feature is stable and relatively finalized, it is time to add it to the documentation. If you are adding any private/public functions, it is best to add docstrings, to aid in reviewing code and also for the API reference.

We use [Numpy style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html>) and [Sphinx](http://www.sphinx-doc.org/en/stable) to document this library. Sphinx, in turn, uses [reStructuredText](http://www.sphinx-doc.org/en/stable/rest.html) as its markup language for adding code.

We use the [Sphinx Autosummary extension](http://www.sphinx-doc.org/en/stable/ext/autosummary.html) to generate API references. You may want to look at `docs/api-*.rst` files to see how they look and where to add new functions, classes or modules. We also use the [myst_nb extension](https://myst-nb.readthedocs.io/en/latest/) to render Jupyter notebooks in the documentation.

To build the documentation, run `sphinx-autobuild` using:

```sh
hatch run docs
```

This will build the documentation and serve it on a local http server which listens for changes and automatically rebuilds.

Documentation from the `main` branch and tagged releases is automatically built and hosted on [readthedocs](https://readthedocs.org/).


## Acknowledgments

This document is based off of the [guidelines from the sparse project](https://github.com/pydata/sparse/blob/master/docs/contributing.rst).
