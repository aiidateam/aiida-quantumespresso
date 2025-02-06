# -*- coding: utf-8 -*-
"""Calcfunction to compile a complete spectrum for each element from multiple powder sample spectra."""
import warnings

from aiida import orm
from aiida.engine import calcfunction
import numpy as np

warnings.warn(
    'This module is deprecated and will be removed soon as part of migrating XAS and XPS workflows to a new repository.'
    '\nThe new repository can be found at: https://github.com/aiidaplugins/aiida-qe-xspec.', FutureWarning
)


@calcfunction
def get_spectra_by_element(elements_list, equivalent_sites_data, **kwargs):
    """Generate a final spectrum for each element from a dictionary of powder spectra inputs.

    Powder spectra to be processed must be passed in using ``kwargs``, in which the keys must
    correspond to the keys of ``equivalent_sites_data``.

    :param elements_list: a List object defining the elements to compile spectra for.
    :param equivalent_sites_data: a Dict object, defining the symmetry properties of the sites associated with each
        powder spectrum in ``kwargs``. Must be in the format used in the ``equivalent_sites_data`` dictionary of
        ``get_xspectra_structures.outputs.output_parameters``
    """

    incoming_spectra_nodes = {key: value for key, value in kwargs.items() if key != 'metadata'}
    elements = elements_list.get_list()
    equivalency_data = equivalent_sites_data.get_dict()
    core_work_chains = {key: value.creator.caller for key, value in incoming_spectra_nodes.items()}

    data_dict = {element: {} for element in elements}
    for key in incoming_spectra_nodes:
        core_work_chain = core_work_chains[key]
        xspectra_out_params = core_work_chain.outputs.parameters_xspectra__xas_0.get_dict()
        energy_zero = xspectra_out_params['energy_zero']
        multiplicity = equivalency_data[key]['multiplicity']
        element = equivalency_data[key]['symbol']

        if 'total_multiplicity' not in data_dict[element]:
            data_dict[element]['total_multiplicity'] = multiplicity
        else:
            data_dict[element]['total_multiplicity'] += multiplicity

        data_dict[element][key] = {
            'spectrum_node': incoming_spectra_nodes[key],
            'element': element,
            'multiplicity': multiplicity,
            'energy_zero': energy_zero
        }

    spectra_by_element = {}
    for element in elements:
        spectra_list = []
        total_multiplicity = data_dict[element].pop('total_multiplicity')
        for key in data_dict[element]:
            spectrum_node = data_dict[element][key]['spectrum_node']
            site_multiplicity = data_dict[element][key]['multiplicity']
            spectrum_x = spectrum_node.get_x()[1]
            spectrum_y = spectrum_node.get_y()[0][1]
            weighted_spectrum = np.column_stack((spectrum_x, (spectrum_y * site_multiplicity) / total_multiplicity))
            spectra_list.append((weighted_spectrum, float(data_dict[element][key]['energy_zero'])))

        # Sort according to Fermi level, then correct to align all spectra to the
        # highest value. Note that this is needed because XSpectra automatically aligns the
        # final spectrum such that the system's Fermi level is at 0 eV.
        spectra_list.sort(key=lambda entry: entry[1])
        highest_level = spectra_list[0][-1]
        energy_zero_corrections = [(entry[0], entry[1] - highest_level) for entry in spectra_list]
        corrected_spectra = [
            np.column_stack((entry[0][:, 0] - entry[1], entry[0][:, 1])) for entry in energy_zero_corrections
        ]

        spectra_by_element[element] = np.column_stack((
            sum(array[:, 0] for array in corrected_spectra) / len(corrected_spectra),
            sum(array[:, 1] for array in corrected_spectra)
        ))

    all_final_spectra = {}
    for element in elements:
        final_spectra = orm.XyData()
        corrected_spectrum = spectra_by_element[element]
        final_spectra_y_labels = [f'{element}_dipole']
        final_spectra_y_units = ['sigma']
        final_spectra_y_arrays = [corrected_spectrum[:, 1]]

        final_spectra_x_label = 'energy'
        final_spectra_x_units = 'eV'
        final_spectra_x_array = corrected_spectrum[:, 0]
        final_spectra.set_x(final_spectra_x_array, final_spectra_x_label, final_spectra_x_units)
        final_spectra.set_y(final_spectra_y_arrays, final_spectra_y_labels, final_spectra_y_units)
        all_final_spectra[f'{element}_xas'] = final_spectra

    return all_final_spectra
