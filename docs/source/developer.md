# Developer Guide

## Quick start

Thanks for contributing! 🙏
Below you can find an overview of the commands to get you up and running quickly, depending on your preferred development tools.

First clone the repository from GitHub and install the package locally in **editable** mode (`-e`):

```
$ git clone https://github.com/aiidateam/aiida-quantumespresso
$ cd aiida-quantumespresso
$ pip install -e .
```

::::{tab-set}

:::{tab-item} Default

The "default" approach to developing is to install the development extras in your current environment:

```console
$ pip install -e .[pre-commit,tests,docs]
```

````{tip}
You can also install the development dependencies using the new [dependency groups](https://packaging.python.org/en/latest/specifications/dependency-groups/), [available from `pip` v25.1](https://ichard26.github.io/blog/2025/04/whats-new-in-pip-25.1/):

```console
$ pip install -e . --group dev
```

In the future the development extras will most likely be removed.

````

**Pre-commit**

To make sure your changes adhere to our formatting/linting preferences, install the pre-commit hooks:

```console
$ pre-commit install
```

They will then run on every `git commit`.
You can also run them on e.g. all files using:

```console
$ pre-commit run -a
```

Drop the `-a` option in case you only want to run on staged files.

**Tests**

You can run all tests in the `tests` directory with `pytest`:

```console
$ pytest
```

Or select the test module:

```console
$ pytest tests/parsers/test_pw.py
```

See the [`pytest` documentation](https://docs.pytest.org/en/stable/how-to/usage.html#specifying-which-tests-to-run) for more information.

**Documentation**

The current documentation build relies on a [Sphinx](https://www.sphinx-doc.org/en/master/index.html)-generated [Makefile](https://www.gnu.org/software/make/manual/make.html).
Build the documentation with:

```console
$ make -C docs html
```

Once complete, you can open the documentation in your browser using:

```console
$ make -C docs view
```

Or clean the documentation files:

```console
$ make -C docs clean
```

See the [documentation](#developer-documentation) section for more information on how we write our documentation.

:::

:::{tab-item} uv

`uv` is a Python package and project manager.
See [the documentation](https://docs.astral.sh/uv/getting-started/installation/) on how to install `uv`.

**Pre-commit**

To make sure your changes adhere to our formatting/linting preferences, install the pre-commit hooks:

```console
$ uvx pre-commit install
```

They will then run on every `git commit`.
You can also run them on e.g. all files using:

```console
$ uvx pre-commit run -a
```

Drop the `-a` option in case you only want to run on staged files.

```{note}
Here we use the [`uvx` command](https://docs.astral.sh/uv/guides/tools/#running-tools) to run the `pre-commit` tool without installing it.
Alternatively you can also install [`pre-commit` as a tool](https://docs.astral.sh/uv/guides/tools/#installing-tools) and omit `uvx`.
```

**Tests**

You can run all tests in the `tests` directory with `pytest`:

```console
$ uv run pytest
```

Or select the test module:

```console
$ uv run pytest tests/parsers/test_pw.py
```

See the [`pytest` documentation](https://docs.pytest.org/en/stable/how-to/usage.html#specifying-which-tests-to-run) for more information.

**Documentation**

The current documentation build relies on a [Sphinx](https://www.sphinx-doc.org/en/master/index.html)-generated [Makefile](https://www.gnu.org/software/make/manual/make.html).
Build the documentation with:

```console
$ uv run make -C docs html
```

Once complete, you can open the documentation in your browser using:

```console
$ uv run make -C docs view
```

Or clean the documentation files:

```console
$ uv run make -C docs clean
```

See the [documentation](#developer-documentation) section for more information on how we write our documentation.

:::

:::{tab-item} Hatch

You can use [Hatch](https://hatch.pypa.io/1.9/install/) to run development tools in isolated environments.

**Pre-commit**

To make sure your changes adhere to our formatting/linting preferences, install the pre-commit hooks:

```console
$ hatch run pre-commit:install
```

They will then run on every `git commit`.
You can also run them on e.g. all files using:

```console
$ hatch run pre-commit:run -a
```

Drop the `-a` option in case you only want to run on staged files.

**Tests**

You can run all tests in the `tests` directory using:

```console
$ hatch test
```

Or select the test module:

```console
$ hatch test tests/parsers/test_pw.py
```

You can also run the tests for a specific Python version with the `-py` option:

```
$ hatch test -py 3.11
```

Or all supported Python `versions` with `--all`:

```
$ hatch test --all
```

See the [Hatch documentation](https://hatch.pypa.io/1.12/tutorials/testing/overview/) for more information.

**Documentation**

The current documentation build relies on a [Sphinx](https://www.sphinx-doc.org/en/master/index.html)-generated [Makefile](https://www.gnu.org/software/make/manual/make.html).
Build the documentation with Hatch in the `docs` environment using:

```console
$ hatch run docs:build
```

This runs `make` under the hood.
Once complete, you can open the documentation in your browser using:

```console
$ hatch run docs:view
```

Or clean the documentation files:

```console
$ hatch run docs:clean
```

See the [documentation](#developer-documentation) section for more information on how we write our documentation.

:::

::::


(developer-documentation)=

## Documentation

We use the [Diátaxis](https://diataxis.fr/) approach for organising the documentation in four sections:

* Tutorials
* How-to's
* Topics (Explanation)
* Reference

All of our documentation should be written in [MyST Markdown](https://myst-parser.readthedocs.io/en/latest).
Below you can find a list of current style guide items:

1. Write **one sentence per line** and otherwise **no manual line wrapping** to make easy to create and review  diffs.
   All standard editors allow for dynamic line wrapping, and the line length is irrelevant for the rendered documentation in, e.g., HTML or PDF format.
1. **File and directory names should be alphanumeric** and all lower-case with underscores as word-separators. Example: `entry_points.rst`
1. **Headers must be set in sentence-case**.
   Example: "Entry points"
1. Separate paragraphs by one empty line, but not more.
1. Use the `-` symbol for itemized lists.

### Notes

#### `run` vs `submit`

Originally discussed in [this issue](https://github.com/aiidateam/aiida-quantumespresso/issues/1127).

In the majority of real use cases, `submit` is the preferred engine command, and the documentation should reflect that.
Hence, we limit using the `run` commando to the quick start documentation, and perhaps showing other use cases such as testing calculation setups or workflows (with caching).
