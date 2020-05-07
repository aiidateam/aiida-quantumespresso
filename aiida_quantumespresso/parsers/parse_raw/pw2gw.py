# -*- coding: utf-8 -*-
"""A collection of function that are used to parse the output of Quantum Espresso pw2gw.

The function that needs to be called from outside is parse_raw_output(). The functions mostly work without aiida
specific functionalities. The parsing will try to convert whatever it can in some dictionary, which by operative
decision doesn't have much structure encoded, [the values are simple ]
"""
from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.parse_raw import convert_qe_time_to_sec
from aiida_quantumespresso.utils.mapping import get_logging_container


def parse_stdout(stdout):
    """Parses the stdout content of a Quantum ESPRESSO `pw2gw.x` calculation.

    :param stdout: the stdout content as a string
    :param input_parameters: dictionary with the input parameters
    :param parser_options: the parser options from the settings input parameter node
    :returns: tuple of two dictionaries, with the parsed data and log messages, respectively
    """
    # Separate the input string into separate lines
    data_lines = stdout.split('\n')

    logs = get_logging_container()

    parsed_data = {}

    for line in data_lines:
        if 'JOB DONE' in line:
            break
    else:
        logs.error.append('ERROR_OUTPUT_STDOUT_INCOMPLETE')

    for count, line in enumerate(data_lines):
        if 'PW2GW' in line and 'WALL' in line:
            try:
                time = line.split('CPU')[1].split('WALL')[0]
                parsed_data['wall_time'] = time
            except Exception:
                logs.warning.append('Error while parsing wall time.')
            try:
                parsed_data['wall_time_seconds'] = convert_qe_time_to_sec(time)
            except ValueError:
                raise QEOutputParsingError('Unable to convert wall_time in seconds.')

    return parsed_data, logs
