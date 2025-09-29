#!/usr/bin/env runaiida
from aiida.engine import run
from aiida.orm import Dict, KpointsData, load_code

# Load the code configured for ``ph.x``. Make sure to replace this string
# with the label of a ``Code`` that you configured in your profile.
code = load_code('ph@localhost')
builder = code.get_builder()

# Replace ``IDENTIFIER_PW_CALCULATION`` with the pk of the completed ``PwCalculation``
builder.parent_folder = load_node(IDENTIFIER_PW_CALCULATION).outputs.remote_folder
builder.parameters = Dict({'INPUTPH': {}})

# Generate a 1x1x1 q-point mesh
qpoints = KpointsData()
qpoints.set_kpoints_mesh([1, 1, 1])
builder.qpoints = qpoints

# Run the calculation on 1 CPU and kill it if it runs longer than 1800 seconds.
# Set ``withmpi`` to ``False`` if ``ph.x`` was compiled without MPI support.
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
