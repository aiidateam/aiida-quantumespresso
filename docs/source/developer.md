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

## Release

Creating a new release consists of the steps outlined below.
In the following, we assume that the latest release was `v3.2.1`.
Furthermore, we assume the following naming conventions for remotes:

* `origin`: remote pointing to `aiidateam/aiida-quantumespresso`
* `fork`: remote pointing to personal fork, e.g. `mbercx/aiida-quantumespresso`

Check how your remotes are configured with `git remote -v`

### Preparing the release branch

#### Deciding the type of release

We use [semantic versioning](https://semver.org/), i.e. version labels have the form `v<major>.<minor>.<patch>`

* Patch release: `v3.2.1` to `v3.2.2`, only bug fixes
* Minor release: `v3.2.1` to `v3.3.0`, bug fixes and new features that **maintain** backwards compatibility
* Major release: `v3.2.1` to `v4.0.0`, bug fixes and new features that **break** backwards compatibility

#### Creating the release branch

We use the [GitHub Flow](https://www.flagship.io/git-branching-strategies/) branching model. In short: new features, bug fixes, etc are added by opening a pull request to `main`, and the `main` branch is tagged after doing a pull request that updates the `CHANGELOG.md`. Doing so will trigger an automated deployment to PyPI.

For most releases, we assume that all the changes in the current `main` branch have to be included in the release. As such, branch off the new release branch directly from `main`:

    git fetch --all --prune
    git checkout origin/main -b release/3.3.0

#### Updating the `CHANGELOG.md`, version and compatibilities

The next step is to update the `CHANGELOG.md` with all the changes made since the last release.
First, update the source code version in the following file by hand:

- `src/aiida_quantumespresso/__init__.py` 

Then, run the `update_changelog.py` script:

```console
python .github/workflows/update_changelog.py
```

This will automatically add:

1. The header for the new release and the sub-headers in which the commits should be sorted.
2. The list of commits since the previous release, with links to them.

Sort the commit items into the correct subsection of the change log.
Ideally, the main changes or new features are also described at the top of the change log message, providing code snippets where it's useful.

Some final checks you should perform:

* If a certain commit introduces or changes functionality, it should be documented.
  If the documentation was not introduced during review (which should always be requested), ask the author of the changes to provide it.

* Also check the [compatibility matrix in `README.md`](https://github.com/aiidateam/aiida-quantumespresso?tab=readme-ov-file#compatibility-matrix) to see if any changes have to be made.

Once you've prepared the release branch locally, commit with the message 'Release `v3.3.0`' and push it to Github. For our major/minor release example:

    git commit -am 'Release `v3.3.0`'
    git push origin release/3.3.0

Your branch is now ready to be released!

### Updating `main`

#### Merge release branch into `main`

Now that the release branch is ready, merge it into `main` via a pull request.
Make sure the remote `release/3.3.0` branch is up to date by pushing the local changes, then go to Github and create a pull request to the `main` branch of the official repository.

After the pull request has been approved, merge the PR using the "Squash and Merge", and make sure the commit is simply named e.g. 'Release `v3.3.0`'.

Once this is complete, fetch all the changes locally, checkout the `main` branch and make sure it's up to date with the remote:

    git fetch --all
    git checkout main
    git pull origin

Next, tag the final release commit:

    git tag -a v3.3.0 -m 'Release `v3.3.0`'

**IMPORTANT**: once you push the tag to GitHub, a workflow will start that automatically publishes a release on PyPI.
Double check that the tag is correct, and that the `main` branch looks in good shape.

If you accidentally tagged the wrong commit, you can delete the local tag using the following command:

    git tag -d v3.3.0

Once you're confident that the tag and `main` branch are in good shape, push both to the remote:

    git push origin main --tags

With the release tag created, the new release is automatically built and published on PyPI via our continuous deployment (CD) GitHub workflow!
