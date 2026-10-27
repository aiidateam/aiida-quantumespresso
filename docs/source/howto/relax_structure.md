(howto-workflows-pw-relax)=

# Relax a structure


The `PwRelaxWorkChain` is designed to perform robust geometry relaxations of atomic structures with Quantum ESPRESSO’s `pw.x`. It combines automated restarts with built‑in error handling.

---

## Minimal example: build and submit

```python
from aiida import orm, load_profile
from aiida_quantumespresso.workflows.pw.relax import PwRelaxWorkChain
from aiida_quantumespresso.common.types import RelaxType
from aiida.engine import submit
from ase.build import bulk

load_profile()

# Load your pw.x code and structure
code = orm.load_code('pw@localhost')
structure = orm.StructureData(ase=bulk('Si', 'diamond', 5.4))

# Get a builder with sensible default parameters from a predefined protocol.
builder = PwRelaxWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure,
    protocol="fast",   # choose from: fast, moderate, precise
    relax_type=RelaxType.POSITIONS_CELL,  # relax both positions and cell
    options={"resources": {"num_machines": 2},
             "max_wallclock_seconds": 7200}
)

# Optionally adjust calculation parameters
builder.base_relax.pw.parameters['SYSTEM']['ecutwfc'] = 60.0
builder.base_relax.pw.parameters['SYSTEM']['ecutrho'] = 300.0

# Submit the work chain
workchain_node = submit(builder)
print(f"Launched {workchain_node.process_label} with PK = {workchain_node.pk}")
```
The available options for `RelaxType`:

* `POSITIONS_CELL`: (default) Optimise both the atomic positions and unit cell.
* `POSITIONS`: Only the atomic positions are relaxed, cell is fixed.
* `SHAPE`: Only the cell shape is optimized at a fixed volume and fixed atomic positions.
* `CELL`: Only the cell is optimized, both shape and volume, while atomic positions are fixed.
* `POSITIONS_SHAPE`: Same as `SHAPE`  but atomic positions are relaxed as well.


---

## Inspecting results

After the work chain finishes, the main outputs are:

* `output_structure`: the relaxed structure.
* `output_parameters`: a dictionary with key results from Quantum ESPRESSO (total energy, Fermi energy, etc.)

```python

print("Work chain finished ok:", workchain_node.is_finished_ok)

# Get the relaxed structure
relaxed_structure = workchain_node.outputs.output_structure
print("Final cell:", relaxed_structure.cell)

# Inspect parsed parameters
params = workchain_node.outputs.output_parameters.get_dict()
print("Total energy:", params.get("energy"))
print("Fermi energy:", params.get("fermi_energy"))
```

You can also inspect the node in the AiiDA database:

```bash
verdi process show <PK>
```
