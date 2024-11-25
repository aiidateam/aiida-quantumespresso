"""A basic parser for the common format of QE."""

import re

from aiida.orm import StructureData as LegacyStructureData
from aiida.plugins import DataFactory
from aiida.common import exceptions

try:
    StructureData = DataFactory('atomistic.structure')
    HAS_ATOMISTIC = True
except exceptions.MissingEntryPointError:
    HAS_ATOMISTIC = False


__all__ = (
    'convert_qe_time_to_sec',
    'convert_qe_to_aiida_structure',
    'convert_qe_to_kpoints',
)


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

    return float(seconds) + float(minutes) * 60.0 + float(hours) * 3600.0 + float(days) * 86400.0


def convert_qe_to_aiida_structure(output_dict, input_structure=None):
    """Convert the dictionary parsed from the Quantum ESPRESSO output into ``StructureData``.
    If we have an ``orm.StructureData`` as input, we return an ``orm.StructureData`` instance,
    otherwise we always return an aiida-atomistic ``StructureData``.
    """

    cell_dict = output_dict['cell']

    # Without an input structure, try to recreate the structure from the output
    if not input_structure:

        if not HAS_ATOMISTIC:
            structure = LegacyStructureData()
            structure.set_cell(cell_dict['lattice_vectors'])

            for kind_name, position in output_dict['atoms']:
                symbol = re.sub(r'\d+', '', kind_name)
                structure.append_atom(position=position, symbols=symbol, name=kind_name)

        else:
            from aiida_atomistic import StructureDataMutable
            structure = StructureDataMutable()
            structure.set_cell(cell_dict['lattice_vectors'])

            for kind_name, position in output_dict['atoms']:
                symbol = re.sub(r'\d+', '', kind_name)
                structure.append_atom(positions=position, symbols=symbol, kinds=kind_name)

            structure = StructureData.from_mutable(structure)

    else:

        if isinstance(input_structure, LegacyStructureData):
            structure = input_structure.clone()
            structure.reset_cell(cell_dict['lattice_vectors'])
            new_pos = [i[1] for i in cell_dict['atoms']]
            structure.reset_sites_positions(new_pos)
        elif HAS_ATOMISTIC:
            if isinstance(input_structure, StructureData):
                structure = input_structure.get_value() # gives the StructureDataMutable instance
                structure.set_cell(cell_dict['lattice_vectors'])
                for site,position in zip(structure.sites,[i[1] for i in cell_dict['atoms']]):
                    site.position = position
            elif isinstance(input_structure, LegacyStructureData):
                structure = input_structure.clone()
                structure.reset_cell(cell_dict['lattice_vectors'])
                new_pos = [i[1] for i in cell_dict['atoms']]
                structure.reset_sites_positions(new_pos)
        else:
            raise ValueError('input_structure is not a valid StructureData or LegacyStructureData instance')
                

    return structure


def convert_qe_to_kpoints(xml_dict, structure):
    """Build the output kpoints from the raw parsed data.

    :param parsed_parameters: the raw parsed data
    :return: a `KpointsData` or None
    """
    from aiida.plugins import DataFactory

    KpointsData = DataFactory('core.array.kpoints')  # noqa: N806

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
