# -*- coding: utf-8 -*-
"""
A basic parser for the common format of QE
"""
from __future__ import absolute_import
from aiida_quantumespresso.parsers import parse_QE_errors, convert_qe_time_to_sec
from aiida_quantumespresso.parsers import QEOutputParsingError, get_parser_info


def parse_qe_simple(filecontent, codename=None):
    """
    Parses the output file of a QE calculation, just checking for basic content
    like JOB DONE, errors with %%%% etc.

    :param filecontent: a string with the output file content
    :param codename: the string printed both in the header and near the walltime.
        If passed, a few more things are parsed (e.g. code version, walltime, ...)
    :return: (successful, out_dict) where successful is a boolean (False is a critical error occurred);
        out_dict is a dictionary with parsed information (e.g. a list of warnings) that could e.g. be
        returned as a Dict by the parser.
    """
    # suppose at the start that the job is successful
    successful = True
    parser_info = get_parser_info(parser_info_template='aiida-quantumespresso parser simple v{}')
    parsed_data = {'warnings': []}
    parsed_data.update(parser_info)

    generic_error_message = "There was an error, please check the 'error_message' key"

    if "JOB DONE" not in filecontent:
        successful = False
        msg = "Computation did not finish properly"
        parsed_data['warnings'].append(msg)

    lines = filecontent.split('\n')

    if codename is not None:
        for count,line in enumerate(lines):

            codestring = "Program {}".format(codename)
            if codestring in line and "starts on" in line:
                parsed_data['code_version'] = line.split(codestring)[1].split("starts on")[0].strip()

            # parse the global file, for informations that are written only once
            if codename in line and 'WALL' in line:
                try:
                    time = line.split('CPU')[1].split('WALL')[0].strip()
                    parsed_data['wall_time'] = time
                except (ValueError, IndexError):
                    parsed_data['warnings'].append('Error while parsing wall time.')
                else:
                    try:
                        parsed_data['wall_time_seconds'] = convert_qe_time_to_sec(time)
                    except ValueError:
                        raise QEOutputParsingError("Unable to convert wall_time in seconds.")

            if '%%%%%%%%%%%%%%' in line:
                if generic_error_message not in parsed_data['warnings']:
                    parsed_data['warnings'].append(generic_error_message)
                if 'error_message' not in parsed_data:
                    parsed_data['error_message'] = []
                successful = False
                # Pass count=0 to start from the top of the file (anyway, it's a short file)
                # pass an empty warnings list because we don't have existing warnings
                # (this is used to avoid duplication of errors)
                messages = parse_QE_errors(lines, count=count,
                                           warnings=parsed_data['error_message'])

                # if it found something, add to log
                if len(messages) > 0:
                    parsed_data['error_message'].extend(messages)

    return successful, parsed_data
