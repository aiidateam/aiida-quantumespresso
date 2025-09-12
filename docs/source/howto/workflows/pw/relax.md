(howto-workflows-pw-relax)=

# `PwRelaxWorkChain`


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
node = submit(builder)
print(f"Launched {node.process_label}<{node.pk}>")
```

---

## Inspecting results

After the work chain finishes, the main outputs are:

* `output_structure`: the relaxed structure.
* `output_parameters`: a dictionary with key results from Quantum ESPRESSO (total energy, Fermi energy, etc.)

```python
from aiida.orm import load_node

wc = load_node(node.pk)
print("Work chain finished ok:", wc.is_finished_ok)

# Get the relaxed structure
if "output_structure" in wc.outputs:
    relaxed_structure = wc.outputs.output_structure
    print("Final cell:", relaxed_structure.cell)

# Inspect parsed parameters
params = wc.outputs.output_parameters.get_dict()
print("Total energy:", params.get("energy"))
print("Fermi energy:", params.get("fermi_energy"))
```

You can also inspect the node in the AiiDA database:

```bash
verdi process show <PK>
```
