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

### Codes

To run a Quantum ESPRESSO (QE) code, it should first be set up in AiiDA.
This can be done easily for multiple codes with the `aiida-quantumespresso` command line interface (CLI).
The command below will set up an AiiDA `Code` for several QE binaries on the computer where AiiDA is running (`localhost`):

```console
$ aiida-quantumespresso setup codes localhost pw.x projwfc.x dos.x
```

:::{important}

The command will look for the executables in your `PATH` using the `which` UNIX command.
This can fail, or you may wish to specify a different executable.
In this case, have a look at the dropdown below.

::::{dropdown} Troubleshooting

You can find out the absolute path to e.g. the first `pw.x` executable in the `PATH` using:

```console
which pw.x
```

If this returns nothing or a Quantum ESPRESSO version you don't want to run, you have a few options:

1. Perhaps you have to set the `PATH`, activate a conda environment, or load a module?
   You can specify the `--prepend-text` as a command-line option, or run in interactive mode with `-i`:

   ```console
   $ aiida-quantumespresso setup codes -i localhost pw.x projwfc.x dos.x
   ```

2. Pass the directory of the correct absolute path as the filepath executable:

   ```console
   $ aiida-quantumespresso setup codes --directory /path/to/bin localhost pw.x projwfc.x dos.x
   ```

::::
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

- the startup file of your shell (`.bashrc`, `.zsh`, ...), if AiiDA is installed system-wide
- the [activators](https://virtualenv.pypa.io/en/latest/user_guide.html#activators) of your virtual environment
- a [startup file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables) for your conda environment

:::{important}
After having added the line to the start up script, make sure to restart the terminal or source the script for the changes to take effect.
:::
