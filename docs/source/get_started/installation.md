---
myst:
  substitutions:
    aiidapseudo: '[`aiida-pseudo`](https://aiida-pseudo.readthedocs.io/en/latest/)'
    pip: '[`pip`](https://pip.pypa.io/en/stable/index.html)'
    PyPI: '[PyPI](https://pypi.org/)'
---

(installation)=

# Installation

The Python package can be installed from the Python Package index {{ PyPI }}:

```console
$ pip install aiida-quantumespresso
```

Once the installation is complete, set up a fresh AiiDA profile if you don't have one yet:

```console
$ verdi presto
```

(installation-setup-code)=

## Basic setup

### `pw.x` code

To run a Quantum ESPRESSO code, it should first be setup in AiiDA.
This can be done from the command line interface (CLI) or the Python application programming interface (API).
In this example, we will setup the `pw.x` code that is installed on the computer where AiiDA is running:

::::{tab-set}

:::{tab-item} CLI

To setup a particular Quantum ESPRESSO code, use the ``verdi`` CLI of ``aiida-core``.

```console
$ verdi code create core.code.installed -n --computer localhost --label pw --default-calc-job-plugin quantumespresso.pw --filepath-executable pw.x
```
:::

:::{tab-item} API

To setup particular Quantum ESPRESSO code using the Python API, run the following code in a Python script with `verdi run` or in the `verdi` shell:

```python
from aiida.orm import InstalledCode

computer = load_computer('localhost')
code = InstalledCode(
label='pw',
computer=computer,
filepath_executable='pw.x',
default_calc_job_plugin='quantumespresso.pw',
).store()
```
:::

::::

:::{important}
Using the commands above, you will set up a code that uses the first `pw.x` binary your `PATH`.
You can find out the absolute path to this binary using the `which` command:

```console
which pw.x
```

If this is not the Quantum ESPRESSO version you want to run, pass the correct absolute path as the filepath executable.
:::

For more detailed information, please refer to the AiiDA documentation [on setting up codes](https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html#how-to-setup-a-code).

(installation-setup-pseudopotentials)=

### Pseudo potentials

Many Quantum ESPRESSO codes require pseudo potentials.
The simplest way of installing these is through the {{ aiidapseudo }} plugin package, which is installed as a dependency of `aiida-quantumespresso`.
At a minimum, at least one pseudo potential family should be installed.
We recommend using the [Standard Solid-State Pseudopotentials (SSSP)](https://www.materialscloud.org/discover/sssp/table/efficiency) v1.3 with the PBEsol functional:

```console
$ aiida-pseudo install sssp -v 1.3 -x PBEsol
```

This is also the default used by our protocols.
For more detailed information on installing other pseudo potential families, please refer to the documentation of {{ aiidapseudo }}.

::::{grid} 1 1 1 1
:gutter: 3

:::{grid-item-card} ✅ All done!
:text-align: center
:shadow: md

The next step is to go to the quick start tutorial to learn the basics of using the package.

+++

```{button-ref} quick_start
:ref-type: doc
:click-parent:
:expand:
:color: primary
:outline:

To the quick start tutorial
```
:::

::::

## Advanced

### CLI tab-completion

To enable tab-completion for the command line interface (CLI), execute the following shell command (depending on the shell):

::::{tab-set}

:::{tab-item} bash
```console
$ eval "$(_AIIDA_QUANTUMESPRESSO_COMPLETE=bash_source aiida-quantumespresso)"
```
:::

:::{tab-item} zsh
```console
$ eval "$(_AIIDA_QUANTUMESPRESSO_COMPLETE=zsh_source aiida-quantumespresso)"
```
:::

:::{tab-item} fish
```console
$ eval (env _AIIDA_QUANTUMESPRESSO_COMPLETE=fish_source aiida-quantumespresso)
```
:::

::::

Place this command in your shell or virtual environment activation script to automatically enable tab completion when opening a new shell or activating an environment.
This file is shell specific, but likely one of the following:

- the startup file of your shell (`.bashrc`, `.zsh`, ...), if aiida is installed system-wide
- the [activators](https://virtualenv.pypa.io/en/latest/user_guide.html#activators) of your virtual environment
- a [startup file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables) for your conda environment

:::{important}
After having added the line to the start up script, make sure to restart the terminal or source the script for the changes to take effect.
:::
