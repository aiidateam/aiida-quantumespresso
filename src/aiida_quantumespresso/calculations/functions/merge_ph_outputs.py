# -*- coding: utf-8 -*-
"""merge data from mulitple ph runs called by one PhBase."""
from aiida import orm
from aiida.engine import calcfunction


@calcfunction
def merge_ph_outputs(**kwargs):
    """Calcfunction to merge outputs from multiple `ph.x` calculations with different q-points."""

    # Get the outputs, sorted by label
    outputs = [el[1] for el in sorted(list(kwargs.items()), key=lambda l: l[0])]

    merged = outputs[-1].get_dict()
    merged.pop('num_q_found', 0)
    total_walltime = 0
    number_of_qpoints = 0
    number_irreps = []

    for index, output in enumerate(outputs):

        output = output.get_dict()
        total_walltime += output.get('wall_time_seconds', 0)
        number_irreps.extend(output['number_of_irr_representations_for_each_q'])
        number_of_qpoints += output.pop('num_q_found', 0)

        for key, value in output.items():

            if key.startswith('dynamical_matrix_') and 'mode_symmetry' in value.keys():
                if '0' in key:
                    merged.pop('dynamical_matrix_0', None)
                    merged[f'dynamical_matrix_{index + 1}'] = value
                else:
                    merged[key] = value

    merged['wall_time_seconds'] = total_walltime
    merged['number_of_irr_representations_for_each_q'] = number_irreps
    merged['number_of_qpoints'] = number_of_qpoints
    merged.pop('wall_time', 0)
    merged = orm.Dict(dict=merged)
    return merged
