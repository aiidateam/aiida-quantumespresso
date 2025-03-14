# -*- coding: utf-8 -*-
"""CalcFunction to merge multiple ``XyData`` nodes of calculated XANES spectra into a new ``XyData`` node."""
import warnings

from aiida.engine import calcfunction
from aiida.orm import XyData

warnings.warn(
    'This module is deprecated and will be removed soon as part of migrating XAS and XPS workflows to a new repository.'
    '\nThe new repository can be found at: https://github.com/aiidaplugins/aiida-qe-xspec.', FutureWarning
)


@calcfunction
def merge_spectra(**kwargs):
    """Compile all calculated spectra into a single ``XyData`` node for easier plotting.

    The keyword arguments must be an arbitrary number of ``XyData`` nodes from
    the `output_spectra` of `XspectraCalculation`s, all other `kwargs` will be discarded at
    runtime.

    Returns a single ``XyData`` node where each set of y values is labelled
    according to the polarisation vector used for the `XspectraCalculation`.
    """
    output_spectra = XyData()
    y_arrays_list = []
    y_units_list = []
    y_labels_list = []

    spectra = [node for label, node in kwargs.items() if isinstance(node, XyData)]

    for spectrum_node in spectra:
        calc_node = spectrum_node.creator
        calc_out_params = calc_node.res
        eps_vector = calc_out_params['xepsilon']

        old_y_component = spectrum_node.get_y()
        if len(old_y_component) == 1:
            y_array = old_y_component[0][1]
            y_units = old_y_component[0][2]
            y_arrays_list.append(y_array)
            y_units_list.append(y_units)
            y_labels_list.append(f'sigma_{eps_vector[0]}_{eps_vector[1]}_{eps_vector[2]}')
        elif len(old_y_component) == 3:
            y_tot = old_y_component[0][1]
            y_tot_units = old_y_component[0][2]
            y_tot_label = f'sigma_tot_{eps_vector[0]}_{eps_vector[1]}_{eps_vector[2]}'
            y_arrays_list.append(y_tot)
            y_units_list.append(y_tot_units)
            y_labels_list.append(y_tot_label)

            y_up = old_y_component[1][1]
            y_up_units = old_y_component[1][2]
            y_up_label = f'sigma_up_{eps_vector[0]}_{eps_vector[1]}_{eps_vector[2]}'
            y_arrays_list.append(y_up)
            y_units_list.append(y_up_units)
            y_labels_list.append(y_up_label)

            y_down = old_y_component[2][1]
            y_down_units = old_y_component[2][2]
            y_down_label = f'sigma_down_{eps_vector[0]}_{eps_vector[1]}_{eps_vector[2]}'
            y_arrays_list.append(y_down)
            y_units_list.append(y_down_units)
            y_labels_list.append(y_down_label)

        x_array = spectrum_node.get_x()[1]
        x_label = spectrum_node.get_x()[0]
        x_units = spectrum_node.get_x()[2]

    output_spectra.set_x(x_array, x_label, x_units)
    output_spectra.set_y(y_arrays_list, y_labels_list, y_units_list)

    return output_spectra
