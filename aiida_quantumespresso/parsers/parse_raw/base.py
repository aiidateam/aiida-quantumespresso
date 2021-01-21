# -*- coding: utf-8 -*-
"""A basic parser for the common format of QE."""

__all__ = ('parse_output_base', 'parse_output_error', 'convert_qe_time_to_sec', 'convert_qe2aiida_structure')


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


def convert_qe2aiida_structure(output_dict, input_structure=None):
    """Receives the dictionary cell parsed from quantum espresso Convert it into an AiiDA structure object."""
    from aiida.plugins import DataFactory

    StructureData = DataFactory('structure')

    cell_dict = output_dict['cell']

    # If I don't have any help, I will set up the cell as it is in QE
    if not input_structure:

        s = StructureData(cell=cell_dict['lattice_vectors'])
        for atom in cell_dict['atoms']:
            s.append_atom(position=tuple(atom[1]), symbols=[atom[0]])

    else:

        s = input_structure.clone()
        s.reset_cell(cell_dict['lattice_vectors'])
        new_pos = [i[1] for i in cell_dict['atoms']]
        s.reset_sites_positions(new_pos)

    return s
