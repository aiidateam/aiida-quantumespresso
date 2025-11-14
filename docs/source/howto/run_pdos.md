(howto-workflows-pdos)=

# Calculate the (projected) DOS

The `PdosWorkChain` is designed to compute the total and projected density of states (DOS and PDOS) for a given structure using Quantum ESPRESSO.
This workflow automates the sequence of calculations required and provides built-in error handling and memory management options.

|                     |                                                               |
|---------------------|---------------------------------------------------------------|
| Workflow class      | {class}`aiida_quantumespresso.workflows.pdos.PdosWorkChain`   |
| Workflow entry point| ``quantumespresso.pdos``                                      |

---

## Overview

Computing DOS and PDOS requires a sequence of four calculations:

1. **SCF calculation** (`pw.x`): Generates the initial self-consistent wavefunction (optional, can be skipped if you provide a parent folder).
2. **NSCF calculation** (`pw.x`): Computes eigenvalues on a denser k-point mesh.
3. **DOS calculation** (`dos.x`): Generates the total density of states from the NSCF results.
4. **PDOS calculation** (`projwfc.x`): Computes the projected density of states by projecting wavefunctions onto atomic orbitals.

The `PdosWorkChain` handles this sequence automatically and provides options for memory management when dealing with large wavefunction files.

---

## Minimal example: build and submit

```python
from aiida import orm, load_profile
from aiida_quantumespresso.workflows.pdos import PdosWorkChain
from aiida.engine import submit
from ase.build import bulk

load_profile()

# Load your codes
pw_code = orm.load_code('pw@localhost')
dos_code = orm.load_code('dos@localhost')
projwfc_code = orm.load_code('projwfc@localhost')

# Create a structure
structure = orm.StructureData(ase=bulk('Si', 'diamond', 5.4))

# Get a builder with sensible default parameters from a predefined protocol
builder = PdosWorkChain.get_builder_from_protocol(
    pw_code=pw_code,
    dos_code=dos_code,
    projwfc_code=projwfc_code,
    structure=structure,
    protocol="balanced",  # choose from: fast, balanced, stringent
    options={
        "resources": {"num_machines": 1},
        "max_wallclock_seconds": 7200
    }
)

# Note that the dos.x and projwfc.x codes are not cpu parallelized so one has to specify the number of mpi processes as 1.
resources = {"num_machines": 1, "num_mpiprocs_per_machine": 1}

builder.projwfc.metadata.options.resources = resources
builder.dos.metadata.options.resources = resources

# Submit the work chain
workchain_node = submit(builder)
print(f"Launched {workchain_node.process_label} with PK = {workchain_node.pk}")
```
For more details on the available protocols and their parameters, see the [protocols topic guide](../topics/protocol).

---

## Skipping the SCF calculation

If you already have a completed SCF calculation and want to reuse the calculated wavefunctions, you can skip the SCF step by not providing the `scf` namespace:

```python
from aiida import orm

# Load a completed SCF calculation
scf_calc = orm.load_node(<PK_OF_SCF_CALCULATION>)

# Build the workflow without the scf namespace
builder = PdosWorkChain.get_builder_from_protocol(
    pw_code=pw_code,
    dos_code=dos_code,
    projwfc_code=projwfc_code,
    structure=structure,
    protocol="moderate"
)

del builder.scf
# Provide the parent folder from the completed SCF
builder.nscf.pw.parent_folder = scf_calc.outputs.remote_folder

builder.projwfc.metadata.options.resources = resources
builder.dos.metadata.options.resources = resources
```

---

## Managing storage memory

The wavefunction files created by the NSCF calculation can become very large (>100 GB), which can cause storage issues when these files are copied to the DOS and PDOS calculations.

### Serial execution with automatic cleanup

Setting `serial_clean` to `True` runs the DOS and PDOS calculations sequentially and automatically cleans up intermediate directories:

```python
builder.serial_clean = True
```

This will:
1. Run the SCF workchain
2. Run the NSCF workchain, then clean the SCF calculation directories
3. Run the DOS calculation, then clean its directory
4. Run the PDOS calculation, then clean its directory

### Clean all working directories after completion

Setting `clean_workdir` to `True` will clean all remaining remote directories after the workflow completes:

```python
builder.clean_workdir = True
```

---

## Inspecting results

After the work chain finishes successfully, the main outputs are:

* `nscf.output_band`: Band structure data from the NSCF calculation
* `nscf.output_parameters`: Parameters from the NSCF calculation 
* `dos.output_dos`: Total density of states data
* `projwfc.projections`: Projections onto atomic orbitals


---


:::{warning}
The `emin`, `emax`, and `deltae` values in the `dos` and `projwfc` inputs must match for the workflow to run correctly.
:::

