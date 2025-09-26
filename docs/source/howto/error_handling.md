(howto-workflows-pw-base)=

# Run with error handling

The `PwBaseWorkChain` is a simple wrapper workflow that runs the `pw.x` code of Quantum ESPRESSO and automatically deals with known errors.


## Basic Example

```python
from ase.build import bulk
from aiida import orm, engine, load_profile
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

load_profile()

code = orm.load_code('pw@localhost')
structure = orm.StructureData(ase=bulk('Si', 'diamond', 5.4))

builder = PwBaseWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure,
    protocol='fast'
)

workchain_node = engine.submit(builder)
```

Once the workflow has completed, you can check the outputs.
For example the `output_parameters` has the energy of the system (in eV):

```python
workchain_node.outputs.output_parameters['energy']
```
