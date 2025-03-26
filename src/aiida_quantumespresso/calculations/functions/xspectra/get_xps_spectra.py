# -*- coding: utf-8 -*-
"""CalcFunction to compute the spectrum from ``XpsWorkchain``."""
import warnings

from aiida import orm
from aiida.engine import calcfunction
import numpy as np

warnings.warn(
    'This module is deprecated and will be removed soon as part of migrating XAS and XPS workflows to a new repository.'
    '\nThe new repository can be found at: https://github.com/aiidaplugins/aiida-qe-xspec.', FutureWarning
)


@calcfunction
def get_spectra_by_element(elements_list, equivalent_sites_data, voight_gamma, voight_sigma, **kwargs):  # pylint: disable=too-many-statements
    """Generate the XPS spectra for each element.

    Calculate the core level shift and binding energy for each element.
    Generate the final spectra using the Voigt profile.

    :param elements_list: a List object defining the list of elements to consider
            when producing spectrum.
    :param equivalent_sites_data: an Dict object containing symmetry data.
    :param voight_gamma: a Float node for the gamma parameter of the voigt profile.
    :param voight_sigma: a Float node for the sigma parameter of the voigt profile.
    :param structure: the StructureData object to be analysed
    :returns: Dict objects for all generated spectra and associated binding energy
            and core level shift.

    """
    from scipy.special import voigt_profile  # pylint: disable=no-name-in-module

    ground_state_node = kwargs.pop('ground_state', None)
    correction_energies = kwargs.pop('correction_energies', orm.Dict()).get_dict()
    incoming_param_nodes = {key: value for key, value in kwargs.items() if key != 'metadata'}
    group_state_energy = None
    if ground_state_node is not None:
        group_state_energy = ground_state_node.get_dict()['energy']
    elements = elements_list.get_list()
    sigma = voight_sigma.value
    gamma = voight_gamma.value
    equivalency_data = equivalent_sites_data.get_dict()

    data_dict = {element: {} for element in elements}
    for key in incoming_param_nodes:
        xspectra_out_params = incoming_param_nodes[key].get_dict()
        multiplicity = equivalency_data[key]['multiplicity']
        element = equivalency_data[key]['symbol']
        total_energy = xspectra_out_params['energy']
        data_dict[element][key] = {'element': element, 'multiplicity': multiplicity, 'total_energy': total_energy}

    result = {}
    core_level_shifts = {}
    binding_energies = {}
    for element in elements:
        spectra_list = []
        for key in data_dict[element]:
            site_multiplicity = data_dict[element][key]['multiplicity']
            spectra_list.append((site_multiplicity, float(data_dict[element][key]['total_energy']), key))
        spectra_list.sort(key=lambda entry: entry[1])
        lowest_total_energy = spectra_list[0][1]
        core_level_shift = [(entry[0], entry[1] - lowest_total_energy, entry[2]) for entry in spectra_list]
        core_level_shifts[element] = core_level_shift
        result[f'{element}_cls'] = orm.Dict(dict={entry[2]: entry[1] for entry in core_level_shift})

        if group_state_energy is not None:
            binding_energy = [(entry[0], entry[1] - group_state_energy + correction_energies[element], entry[2])
                              for entry in spectra_list]
            binding_energies[element] = binding_energy
            result[f'{element}_be'] = orm.Dict(dict={entry[2]: entry[1] for entry in binding_energy})

    fwhm_voight = gamma / 2 + np.sqrt(gamma**2 / 4 + sigma**2)

    def spectra_broadening(points, label='cls_spectra'):
        """Broadening base on the binding energy."""
        result_spectra = {}
        for element in elements:
            final_spectra_y_arrays = []
            final_spectra_y_labels = []
            final_spectra_y_units = []

            total_multiplicity = sum(i[0] for i in points[element])

            final_spectra = orm.XyData()
            max_core_level_shift = points[element][-1][1]
            min_core_level_shift = points[element][0][1]
            # Energy range for the Broadening function
            x_energy_range = np.linspace(
                min_core_level_shift - fwhm_voight - 1.5, max_core_level_shift + fwhm_voight + 1.5, 500
            )

            for atoms, index in zip(points[element], range(len(points[element]))):
                # Weight for the spectra of every atom
                intensity = atoms[0]
                relative_peak_position = atoms[1]
                final_spectra_y_labels.append(f'{element}{index}_xps')
                final_spectra_y_units.append('sigma')
                final_spectra_y_arrays.append(
                    intensity * voigt_profile(x_energy_range - relative_peak_position, sigma, gamma) /
                    total_multiplicity
                )

            final_spectra_y_labels.append(f'{element}_total_xps')
            final_spectra_y_units.append('sigma')
            final_spectra_y_arrays.append(sum(final_spectra_y_arrays))

            final_spectra_x_label = 'energy'
            final_spectra_x_units = 'eV'
            final_spectra_x_array = x_energy_range
            final_spectra.set_x(final_spectra_x_array, final_spectra_x_label, final_spectra_x_units)
            final_spectra.set_y(final_spectra_y_arrays, final_spectra_y_labels, final_spectra_y_units)
            result_spectra[f'{element}_{label}'] = final_spectra
        return result_spectra

    result.update(spectra_broadening(core_level_shifts))
    if ground_state_node is not None:
        result.update(spectra_broadening(binding_energies, label='be_spectra'))
    return result
