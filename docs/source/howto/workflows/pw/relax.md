(howto-workflows-pw-relax)=

# `PwRelaxWorkChain`

As with all work chains, the best way to get started with running the `PwRelaxWorkChain` is *via* the `get_builder_from_protocol()` method.

First, load the AiiDA orm and activate your profile:

```python
from aiida import orm, load_profile

load_profile()
```

Next, load the `PwRelaxWorkChain` and the structure/code you want to run:

```python
from aiida_quantumespresso.workflows.pw.relax import PwRelaxWorkChain

code = orm.load_code('qe-7.2-pw@localhost')
structure = orm.load_node(1)  # FCC Silicon
```

Now you can obtain a fully populated builder using the `get_builder_from_protocol()` method:

```python
builder = PwRelaxWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure
)
```

You can check its contents using:

```python
builder
```

:::{important}
The `get_builder_from_protocol()` method provides you with a fully populated builder with sensible defaults for most structures.
However, it does *not* guarantee convergence for all structures, so it is important to double-check e.g. k-points meshes and plane wave cutoffs for important results.
:::

Now that you have a fully-populated builder that is ready to run, load the `submit` function from the AiiDA engine and launch the process!

```python
from aiida.engine import submit

submit(builder)
```

## Defining the `RelaxType`

By default, the `get_builder_from_protocol()` method will set inputs that instruct Quantum ESPRESSO to optimise both the geometry and cell, i.e.:

* `CONTROL.calculation = vc-relax`
* `CELL.conv_thr = all`

You can of course change these inputs by adapting the `builder` afterwards, for example:

```python
builder = PwRelaxWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure
)

parameters = builder.base.pw.parameters.get_dict()
parameters['CONTROL']['calculation'] = 'relax'
parameters.pop('CELL')
builder.base.pw.parameters = orm.Dict(parameters)
```

```python
from aiida.engine import submit

submit(builder)
```

However, it is somewhat tedious to have to adapt several parameters, and this requires knowledge about the correct combinations of inputs in Quantum ESPRESSO.
To make it easier to specify what degrees of freedom to allow during the geometry optimization, you can also use the `RelaxType` enum:

```python
from aiida_quantumespresso.common.types import RelaxType

builder = PwRelaxWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure,
    relax_type=RelaxType.POSITIONS
)
```

Check the contents of the `builder`!
You'll see that the same inputs as above are now set, so only the _positions_ of the atoms are optimized.
Here is the list of settings for `RelaxType` that are currently supported:

* `POSITIONS`; Only the atomic positions are relaxed, cell is fixed.
* `SHAPE`: Only the cell shape is optimized at a fixed volume and fixed atomic positions.
* `CELL`: Only the cell is optimized, both shape and volume, while atomic positions are fixed.
* `POSITIONS_SHAPE`: Same as `SHAPE`  but atomic positions are relaxed as well.
* `POSITIONS_CELL`: Same as `CELL`  but atomic positions are relaxed as well.

As was noted above, the default setting here is to optimized both the positions and the unit cell, i.e. `RelaxType.POSITIONS_CELL`.
