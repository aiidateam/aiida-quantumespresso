(howto-workflows-pw-relax)=

# `PwRelaxWorkChain`

First, load the AiiDA orm and activate your profile:

```python
from aiida import orm, load_profile

load_profile()
```

Import the `PwRelaxWorkChain` and use the protocols to run a `structure` with the `code` you have set up:

```python
from ase.build import bulk
from aiida_quantumespresso.workflows.pw.relax import PwRelaxWorkChain

code = orm.load_code('pw@localhost')
structure = orm.StructureData(ase=bulk('Si', 'diamond', 5.4))

builder = PwRelaxWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure,
    protocol='fast'
)
```

The calculation can then be executed using for example the `run` function of the engine:

```python
from aiida.engine import run

results = run(builder)
```

The `results` will contain the output nodes, for example the `output_structure`:

```python
results['output_structure'].get_pymatgen()
```

## Define the `RelaxType`

By default, the `get_builder_from_protocol()` method will set inputs that instruct Quantum ESPRESSO to optimise both the geometry and cell, i.e.:

* `CONTROL.calculation = vc-relax`
* `CELL.conv_thr = all`

In case you want to _only_ optimise the positions, you can specify this using the `relax_type` input:


```python
from aiida_quantumespresso.common.types import RelaxType

builder = PwRelaxWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure,
    relax_type=RelaxType.POSITIONS
)
```

Options for `RelaxType`:

* `POSITIONS_CELL`: (default) Optimise both the atomic positions and unit cell.
* `POSITIONS`; Only the atomic positions are relaxed, cell is fixed.
* `SHAPE`: Only the cell shape is optimized at a fixed volume and fixed atomic positions.
* `CELL`: Only the cell is optimized, both shape and volume, while atomic positions are fixed.
* `POSITIONS_SHAPE`: Same as `SHAPE`  but atomic positions are relaxed as well.
