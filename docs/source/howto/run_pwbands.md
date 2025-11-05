(howto-workflows-pw-bands)=

# Calculate a band structure

The `PwBandsWorkChain` is designed to compute electronic band structures for a given structure using Quantum ESPRESSO's `pw.x`.
It automates the complete workflow, starting from a previously relaxed structure then it performs an SCF calculation, and bands calculation along high-symmetry k-points paths.

|                     |                                                               |
|---------------------|---------------------------------------------------------------|
| Workchain class     | {class}`aiida_quantumespresso.workflows.pw.bands.PwBandsWorkChain` |
| Workchain entry point | ``quantumespresso.pw.bands``                                |

---

## Minimal example: build and submit

```python
from aiida import orm, load_profile
from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain
from aiida.engine import submit
from ase.build import bulk

load_profile()

# Load your pw.x code and structure
code = orm.load_code('pw@localhost')
structure = orm.StructureData(ase=bulk('Si', 'diamond', 5.4))

# Get a builder with sensible default parameters from a predefined protocol
builder = PwBandsWorkChain.get_builder_from_protocol(
    code=code,
    structure=structure,
    protocol="balanced",  # choose from: fast, balanced, stringent
    options={
        "account": "your_account", # Change to your account if needed by your HPC provider. Otherwise, remove this line.
        "queue_name": "debug",
        "resources": {"num_machines": 1},
        "max_wallclock_seconds": 1800,
    }
)

# Submit the work chain
workchain_node = submit(builder)
print(f"Launched {workchain_node.process_label} with PK = {workchain_node.pk}")
```

The workchain will automatically:
1. Use SeeKpath to find the primitive cell and generate a high-symmetry k-points path.
2. Run an SCF calculation to obtain the ground state.
3. Run a bands calculation along the high-symmetry k-points path.

---

## Customizing the workflow

### Specifying custom k-points list

You can provide your own k-points instead of using SeeKpath:

```python
from aiida import orm

# Create a custom KpointsData. For example, a list of high-symmetry points:
kpoints = orm.KpointsData()
kpoints.set_kpoints([
    [0.0, 0.0, 0.0],  # Gamma
    [0.5, 0.0, 0.5],  # X
    [0.5, 0.25, 0.75],  # W
    # ... more k-points
])

builder.bands_kpoints = kpoints
del builder.bands_kpoints_distance # Remove the distance input if previously set
```
In that case, you need to remove the `bands_kpoints_distance` input that was set by default by the protocol.

### Adjusting the number of bands

Control the number of bands in the bands calculation using `nbands_factor`:

```python
# Number of bands will be: max(nbnd_scf, 0.5 * nelectrons * factor, 0.5 * nelectrons + 4)
builder.nbands_factor = orm.Float(1.5)
```

Alternatively, set the number of bands directly:

```python
builder.bands.pw.parameters['SYSTEM']['nbnd'] = 50
```

**Note:** You cannot specify both `nbands_factor` and `bands.pw.parameters.SYSTEM.nbnd`.


---

## Inspecting results

After the work chain finishes, the main outputs are:

* `band_structure`: the computed electronic band structure along the k-points path
* `scf_parameters`: output parameters from the SCF calculation
* `band_parameters`: output parameters from the bands calculation
* `primitive_structure`: the normalized primitive structure (if SeeKpath was used)
* `seekpath_parameters`: parameters from the SeeKpath analysis (if SeeKpath was used)


### Visualizing the band structure

You can visualize the band structure using common plotting libraries. For example here, a simple script that shows the general structure of the bands:

```python
import matplotlib.pyplot as plt
import numpy as np

band_structure = workchain_node.outputs.band_structure
seekpath_params = workchain_node.outputs.seekpath_parameters.get_dict()
explicit_kpoints = seekpath_params['explicit_kpoints_linearcoord']
bands = band_structure.get_bands()

for i in range(bands.shape[1]):
     plt.plot(explicit_kpoints, bands[:, i], color='blue')
```

### Protocol details

The available protocols (See the [discussion](#protocols) for more details ) (`fast`, `balanced`, `stringent`) differ in:
- Plane-wave cutoff energies
- k-points density for SCF
- k-points spacing along the bands path
- Convergence thresholds

Choose `fast` for quick tests, `balanced` for production calculations, and `stringent` for high-accuracy results.
