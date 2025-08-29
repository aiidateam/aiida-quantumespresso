---
jupyter:
  jupytext:
    cell_metadata_filter: -all
    formats: .jupytext-sync-ipynb//ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.17.2
  kernelspec:
    display_name: QE
    language: python
    name: python3
---

(howto-workflows-pw-base)=

# `PwBaseWorkChain`

The `PwBaseWorkChain` is a simple wrapper workflow that runs the `pw.x` code of Quantum ESPRESSO and automatically deals with known errors.
To run the code in this how-to, start by loading the AiiDA ORM module and loading the profile:

```python
from aiida import orm, load_profile

load_profile()
```

## Protocols

The fastest way to get started running a `PwBaseWorkChain` is using the protocols:

```python
from ase.build import bulk
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

code = orm.load_code('pw@localhost')
structure = orm.StructureData(ase=bulk('Si', 'diamond', 5.4))

builder = PwBaseWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure,
    protocol='fast'
)
```

The calculation can then be run using for example the `run` function of the engine:

```python
from aiida import engine

workchain_node = engine.submit(builder)
```

With the `fast` protocol, the `pw.x` calculation should be quick.
Once the workflow has completed, you can check the outputs.
For example the `output_parameters` has the energy of the system (in eV):

```python
workchain_node.outputs.output_parameters['energy']
```

```python

```
