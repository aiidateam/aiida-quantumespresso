---
myst:
  substitutions:
    SSSP: Standard Solid-State Pseudopotentials (SSSP)
    aiida-pseudo: '`aiida-pseudo`'
    pip: '[`pip`](https://pip.pypa.io/en/stable/index.html)'
    PyPI: '[PyPI](https://pypi.org/)'
---

# Get started

(installation-requirements)=

## Requirements

To work with `aiida-quantumespresso`, you should have:

- installed `aiida-core`
- configured an AiiDA profile.

Please refer to the [documentation](https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/get_started.html) of `aiida-core` for detailed instructions.

(installation-installation)=

## Installation

The Python package can be installed from the Python Package index {{ PyPI }} or directly from the source:

::::{tab-set}

:::{tab-item} PyPI
The recommended method of installation is to use the Python package manager {{ pip }}:

```console
$ pip install aiida-quantumespresso
```

This will install the latest stable version that was released to {{ PyPI }}.
:::

:::{tab-item} Source
To install the package from source, first clone the repository and then install using {{ pip }}:

```console
$ git clone https://github.com/aiidateam/aiida-quantumespresso
$ pip install -e aiida-quantumespresso
```

The ``-e`` flag will install the package in editable mode, meaning that changes to the source code will be automatically picked up.
:::

::::

(installation-configuration)=

## Configuration

To enable tab-completion for the command line interface, execute the following shell command (depending on the shell):

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

(installation-setup)=

## Setup

(installation-setup-computer)=

### Computer

To run Quantum ESPRESSO calculations on a compute resource, the computer should first be set up in AiiDA.
This can be done from the command line interface (CLI) or the Python application programming interface (API).
In this example, we will set up the `localhost`, the computer where AiiDA itself is running:

::::{tab-set}

:::{tab-item} CLI

To set up a computer, use the ``verdi`` CLI of ``aiida-core``.

```console
$ verdi computer setup -n -L localhost -H localhost -T core.local -S core.direct -w $AIIDA_PATH/workdir
```

After creating the localhost computer, configure the `core.local` transport using:

```console
$ verdi computer configure core.local localhost -n --safe-interval 0
```

Verify that the computer was properly setup by running:

```console
$ verdi computer test localhost
```
:::

:::{tab-item} API

To setup a computer using the Python API, run the following code in a Python script with `verdi run` or in the `verdi` shell:

```python
from aiida.orm import Computer
from pathlib import Path

computer = Computer(
    label='localhost',
    hostname='localhost',
    transport_type='core.local',
    scheduler_type='core.direct',
    workdir=Path('$AIIDA_PATH/workdir').absolute().as_posix()
).store()
computer.configure()
```
:::
::::

For more detailed information, please refer to the documentation [on setting up compute resources](https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html#how-to-set-up-a-computer).

(installation-setup-code)=

### Code

To run a Quantum ESPRESSO code, it should first be setup in AiiDA.
This can be done from the command line interface (CLI) or the Python application programming interface (API).
In this example, we will setup the `pw.x` code that is installed on the computer where AiiDA is running:

::::{tab-set}

:::{tab-item} CLI

To setup a particular Quantum ESPRESSO code, use the ``verdi`` CLI of ``aiida-core``.

```console
$ verdi code create core.code.installed -n --computer localhost --label pw --default-calc-job-plugin quantumespresso.pw --filepath-executable /path/to/pw.x
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
filepath_executable='/path/to/pw.x',
default_calc_job_plugin='quantumespresso.pw',
).store()
```
:::

::::

:::{important}
Make sure to replace `/path/to/pw.x` with the actual absolute path to the `pw.x` binary.
If the binary is in your `PATH`, you can find the first `pw.x` binary in your `PATH` directories using the `which` command:

```console
which pw.x
```
:::

For more detailed information, please refer to the documentation [on setting up codes](https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html#how-to-setup-a-code).

(installation-setup-pseudopotentials)=

### Pseudo potentials

Many Quantum ESPRESSO codes require pseudo potentials.
The simplest way of installing these is through the `aiida-pseudo` plugin package.
This should come as a dependency of `aiida-quantumespresso` and should already be installed.
If this is not the case, it can be installed using:

```console
$ pip install aiida-pseudo
```

At a minimum, at least one pseudo potential family should be installed.
We recommend using the [SSSP] with the PBEsol functional:

```console
$ aiida-pseudo install sssp -x PBEsol
```

For more detailed information on installing other pseudo potential families, please refer to the documentation of [aiida-pseudo].
