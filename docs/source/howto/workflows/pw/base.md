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

results = engine.run(builder)
```

With the `fast` protocol, this calculation should be quick.
The `results` variable contains all the output nodes, for example the `output_parmaters` containing the energy of the system (in eV):

```python
results['output_parameters'].get_dict()['energy']
```
