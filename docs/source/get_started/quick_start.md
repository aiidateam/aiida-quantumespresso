---
jupyter:
  jupytext:
    formats: .jupytext-sync-ipynb//ipynb,md
    main_language: python
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.17.3
---

(quick-start)=

# Quick start

In this tutorial we'll show the basics of running Quantum ESPRESSO calculations with AiiDA.
Start by loading the ORM module and AiiDA profile:

```python
from aiida import orm, load_profile

load_profile()
```

## Running your first `pw.x` calculation

Let's first create a silicon structure - for example with [ASE](https://ase-lib.org/) - and store it in a `StructureData` node.

```python
from ase.build import bulk

structure = orm.StructureData(ase=bulk('Si', 'fcc', 5.43))
```

Next, we'll load the `pw.x` code, which should have been set up in [basic installation instructions](#installation-setup-code).

```python
code = orm.load_code('pw@localhost')
```

The AiiDA Quantum ESPRESSO plugin comes with a set of tested protocols, which define default input parameters for different levels of precision.
You can obtain pre-populated `builder` for most of the supported workflows calling the `get_builder_from_protocol` method:

```python
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

builder = PwBaseWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure,
    protocol='fast'
)
builder
```

In short, the [Process Builder](https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/processes/usage.html#process-builder) is an AiiDA tool that allows you to set up the inputs for a process.
Using the protocols, you get a fully populated one, that can be run immediately by the AiiDA engine:

```python
from aiida import engine

results = engine.run(builder)
```

Because we selected the `fast` protocol, this calculation should finish quickly.
The `results` contain the output nodes of the process.
For example the calculated energy of the system in eV:

```python
results['output_parameters']['energy']
```

## Setting inputs

In many cases you will want to adapt the inputs that are provided by the protocol.
One way to do this is by providing an `overrides` dictionary:

```python
overrides = {
    'pw': {
        'parameters': {
            'SYSTEM': {'nbnd': 10}
        }
    }
}

builder = PwBaseWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure,
    protocol='fast',
    overrides=overrides
)
builder
```

Alternatively, you can also adapt the inputs on the `builder` afterwards:

```python
builder.pw.parameters['SYSTEM']['nbnd'] = 15
builder

```

## Submitting and inspecting processes

So far, we've used the `engine.run()` function to run the `PwBaseWorkChain`.
This is very useful for short running processes in tutorials such as this one, or while developing a workflow.
However, if the `run` command is interrupted, so is the `pw.x` calculation we are running.

Instead, we typically want to `submit` the process to the AiiDA _daemon_, which will take care of its execution.

```python
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

builder = PwBaseWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure,
    protocol='fast'
)
workchain_node = engine.submit(builder)
```

:::{important}

The `submit()` function returns the work chain _node_, not the results!
This makes sense: after submitting the calculation, it hasn't been run yet, so there are no results to return.

:::

The AiiDA daemon is a process that can run in the background and manage your calculations and workflows.
Let's see if it's running:

:::{margin}
**Note**: Lines starting with `!` aren’t Python.
This is an [IPython/Jupyter shell escape](https://ipython.readthedocs.io/en/stable/interactive/reference.html#system-shell-access) — the command after `!` is run in your system shell (e.g. `bash` or `zsh`), not in Python itself.
Alternatively, you can open a terminal and run these commands there (minus the `!`).
:::

```
!verdi daemon status
```

:::{admonition} The `verdi` CLI
AiiDA has a command line interface (CLI): `verdi`.
You can read more about it in the [`aiida-core` documentation](https://aiida.readthedocs.io/projects/aiida-core/en/stable/topics/cli.html).
:::

If it's not, start it:

```
!verdi daemon start
```

Now, let's have a look at the list of processes we have run:

```
!verdi process list -a -p1
```

By default, `verdi process list` will only list "active" processes.
Using the `-a/--all` option, we can see inactive ones too.
Adding `-p/--past-days 1` only shows processes run in the last day, which is a good trick to only list recent processes.

If all is well, you should see that our `PwBaseWorkChain` and the `PwCalculation` it calls are `Finished` with exit status `[0]`.
Zero here means the process finished _successfully_, whereas non-zero numbers typically indicate something went wrong.
You can also inspect the `workchain_node` in the Python API:

```python
workchain_node.is_finished_ok
```

:::{tip}
Restarted your kernel and lost your variables?
Don’t worry — your data is still in the AiiDA database.
You can load a node back into a variable with its UUID or PK:

```python
workchain_node = orm.load_node(<PK>)  # Replace <PK> here!
```

This gives you the same `Node` object again.
:::

To wrap up, let's look at the same result as we did before:

```python
workchain_node.outputs.output_parameters['energy']
```

Hopefully that result is the same as the one above for a reasonable number of digits!
The `output_parameters` are only one of the `outputs` of the `PwBaseWorkChain`:

```
list(workchain_node.outputs)
```

You can also inspect the in- and outputs of the process with the `verdi` CLI:

:::{margin}
Replace `<PK>` here!
:::

```
!verdi process show <PK>
```
