# -*- coding: utf-8 -*-
"""A basic parser for the common format of QE."""
import re

from aiida.orm import StructureData

__all__ = ('convert_qe_time_to_sec', 'convert_qe_to_aiida_structure', 'convert_qe_to_kpoints')


def convert_qe_time_to_sec(timestr):
    """Given the walltime string of Quantum Espresso, converts it in a number of seconds (float)."""
    rest = timestr.strip()

    if 'd' in rest:
        days, rest = rest.split('d')
    else:
        days = '0'

    if 'h' in rest:
        hours, rest = rest.split('h')
    else:
        hours = '0'

    if 'm' in rest:
        minutes, rest = rest.split('m')
    else:
        minutes = '0'

    if 's' in rest:
        seconds, rest = rest.split('s')
    else:
        seconds = '0.'

    if rest.strip():
        raise ValueError(f"Something remained at the end of the string '{timestr}': '{rest}'")

    num_seconds = (float(seconds) + float(minutes) * 60. + float(hours) * 3600. + float(days) * 86400.)

    return num_seconds


def convert_qe_to_aiida_structure(output_dict, input_structure=None):
    """Convert the dictionary parsed from the Quantum ESPRESSO output into ``StructureData``."""

    cell_dict = output_dict['cell']

    # Without an input structure, try to recreate the structure from the output
    if not input_structure:

        structure = StructureData(cell=cell_dict['lattice_vectors'])

        for kind_name, position in output_dict['atoms']:
            symbol = re.sub(r'\d+', '', kind_name)
            structure.append_atom(position=position, symbols=symbol, name=kind_name)

    else:

        structure = input_structure.clone()
        structure.reset_cell(cell_dict['lattice_vectors'])
        new_pos = [i[1] for i in cell_dict['atoms']]
        structure.reset_sites_positions(new_pos)

    return structure


def convert_qe_to_kpoints(xml_dict, structure):
    """Build the output kpoints from the raw parsed data.

    :param parsed_parameters: the raw parsed data
    :return: a `KpointsData` or None
    """
    from aiida.plugins import DataFactory

    KpointsData = DataFactory('core.array.kpoints')

    k_points_list = xml_dict.get('k_points', None)
    k_points_units = xml_dict.get('k_points_units', None)
    k_points_weights_list = xml_dict.get('k_points_weights', None)

    if k_points_list is None or k_points_weights_list is None:
        return None

    if k_points_units != '1 / angstrom':
        raise ValueError('k-points are not expressed in reciprocal cartesian coordinates')

    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    kpoints.set_kpoints(k_points_list, cartesian=True, weights=k_points_weights_list)

    return kpoints
