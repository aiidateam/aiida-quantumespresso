# -*- coding: utf-8 -*-
"""merge data from mulitple ph runs called by one PhBase."""
from aiida import orm
from aiida.engine import calcfunction


@calcfunction
def merge_ph_outputs(**kwargs):
    """Calcfunction to merge outputs from multiple `ph.x` calculations with different q-points."""

    # Get the outputs, sorted by label
    outputs = [el[1].get_dict() for el in sorted(list(kwargs.items()), key=lambda l: l[0])]

    merged = {}

    total_walltime = 0
    number_of_qpoints = 0
    number_irreps = []

    for output in outputs:

        num_irreps_per_q = output.pop('number_of_irr_representations_for_each_q', [])
        number_of_qpoints += len(num_irreps_per_q)
        number_irreps.extend(num_irreps_per_q)
        total_walltime += output.pop('wall_time_seconds', 0)

        for key, value in output.items():
            # Skip uncomplete dynamical matrix results
            if 'dynamical_matrix_' in key and 'mode_symmetry' not in value:
                continue
            merged[key] = value

    merged['wall_time_seconds'] = total_walltime
    merged['number_of_irr_representations_for_each_q'] = number_irreps
    merged['number_of_qpoints'] = number_of_qpoints

    return orm.Dict(dict=merged)
