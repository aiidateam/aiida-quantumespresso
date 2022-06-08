# -*- coding: utf-8 -*-
"""A basic parser for the common format of QE."""
import re

from aiida.orm.nodes.data.structure import Kind, Site
from aiida.plugins import DataFactory

StructureData = DataFactory('core.structure')

__all__ = ('parse_output_base', 'parse_output_error', 'convert_qe_time_to_sec', 'convert_qe_to_aiida_structure')


def parse_output_base(filecontent, codename=None, message_map=None):
    """Parses the output file of a QE calculation, just checking for basic content like JOB DONE, errors with %%%% etc.

    :param filecontent: a string with the output file content
    :param codename: the string printed both in the header and near the walltime.
        If passed, a few more things are parsed (e.g. code version, walltime, ...)
    :returns: tuple of two dictionaries, with the parsed data and log messages, respectively
    """
    from aiida_quantumespresso.utils.mapping import get_logging_container

    keys = ['error', 'warning']

    if message_map is not None and (not isinstance(message_map, dict) or any(key not in message_map for key in keys)):
        raise RuntimeError(f'invalid format `message_map`: should be dictionary with two keys {keys}')

    logs = get_logging_container()
    parsed_data = {}

    lines = filecontent if isinstance(filecontent, list) else filecontent.split('\n')

    for line in lines:
        if 'JOB DONE' in line:
            break
    else:
        logs.error.append('ERROR_OUTPUT_STDOUT_INCOMPLETE')

    if codename is not None:

        codestring = f'Program {codename}'

        for line_number, line in enumerate(lines):

            if codestring in line and 'starts on' in line:
                parsed_data['code_version'] = line.split(codestring)[1].split('starts on')[0].strip()

            # Parse the walltime
            if codename in line and 'WALL' in line:
                try:
                    time = line.split('CPU')[1].split('WALL')[0].strip()
                    parsed_data['wall_time'] = time
                except (ValueError, IndexError):
                    logs.warnings.append('ERROR_PARSING_WALLTIME')
                else:
                    try:
                        parsed_data['wall_time_seconds'] = convert_qe_time_to_sec(time)
                    except ValueError:
                        logs.warnings.append('ERROR_CONVERTING_WALLTIME_TO_SECONDS')

            # Parse an error message with optional mapping of the message
            if '%%%%%%%%%%%%%%' in line:
                parse_output_error(lines, line_number, logs, message_map)

    return parsed_data, logs


def parse_output_error(lines, line_number_start, logs, message_map=None):
    """Parse a Quantum ESPRESSO error message which appears between two lines marked by ``%%%%%%%%``)

    :param lines: a list of strings gotten by splitting the standard output content on newlines
    :param line_number_start: the line at which we identified some ``%%%%%%%%``
    :param logs: a logging container from `aiida_quantumespresso.utils.mapping.get_logging_container`
    """

    def map_message(message, message_map, logs):

        # Match any known error and warning messages
        for marker, message in message_map['error'].items():
            if marker in line:
                if message is None:
                    message = line
                logs.error.append(message)

        for marker, message in message_map['warning'].items():
            if marker in line:
                if message is None:
                    message = line
                logs.warning.append(message)

    # First determine the line that closes the error block which is also marked by ``%%%%%%%`` in the line
    for line_number, line in enumerate(lines[line_number_start + 1:]):
        if '%%%%%%%%%%%%' in line:
            line_number_end = line_number
            break
    else:
        return

    # Get the set of unique lines between the error indicators and pass them through the message map, or if not provided
    # simply append the message to the `error` list of the logs container
    for message in set(lines[line_number_start:line_number_end]):
        if message_map is not None:
            map_message(message, message_map, logs)
        else:
            logs.error(message)

    return


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
