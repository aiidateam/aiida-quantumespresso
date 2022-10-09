#!/usr/bin/env runaiida
# -*- coding: utf-8 -*-
from aiida.engine import run
from aiida.orm import Dict, KpointsData, StructureData, load_code, load_group

# Load the code configured for ``pw.x``. Make sure to replace this string
# with the label that you used in the :ref:`code setup <installation:setup:code>`_.
code = load_code('pw@localhost')
builder = code.get_builder()

# Create a silicon fcc crystal
from ase.build import bulk

structure = StructureData(ase=bulk('Si', 'fcc', 5.43))
builder.structure = structure

# Load the pseudopotential family.
pseudo_family = load_group('SSSP/1.1/PBE/efficiency')
builder.pseudos = pseudo_family.get_pseudos(structure=structure)

# Request the recommended wavefunction and charge density cutoffs
# for the given structure and energy units.
cutoff_wfc, cutoff_rho = pseudo_family.get_recommended_cutoffs(
    structure=structure,
    unit='Ry'
)

parameters = Dict({
    'CONTROL': {
        'calculation': 'scf'
    },
    'SYSTEM': {
        'ecutwfc': cutoff_wfc,
        'ecutrho': cutoff_rho,
    }
})
builder.parameters = parameters

# Generate a 2x2x2 Monkhorst-Pack mesh
kpoints = KpointsData()
kpoints.set_kpoints_mesh([2, 2, 2])
builder.kpoints = kpoints

# Run the calculation on 1 CPU and kill it if it runs longer than 1800 seconds.
# Set ``withmpi`` to ``False`` if ``pw.x`` was compiled without MPI support.
builder.metadata.options = {
    'resources': {
        'num_machines': 1,
    },
    'max_wallclock_seconds': 1800,
    'withmpi': False,
}

results, node = run.get_node(builder)
print(f'Calculation: {node.process_class}<{node.pk}> {node.process_state.value} [{node.exit_status}]')
print(f'Results: {results}')
assert node.is_finished_ok, f'{node} failed: [{node.exit_status}] {node.exit_message}'
