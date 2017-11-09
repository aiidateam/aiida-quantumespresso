# -*- coding: utf-8 -*-
"""
A basic parser for the common format of QE
"""
from aiida_quantumespresso.parsers import parse_QE_errors

def parse_qe_simple(filecontent):
    """
    Parses the output file of a QE calculation, just checking for basic content
    like JOB DONE, errors with %%%% etc.

    :param filecontent: a string with the output file content
    :return: (successful, out_dict) where successful is a boolean (False is a critical error occurred);
        out_dict is a dictionary with parsed information (e.g. a list of warnings) that could e.g. be
        returned as a ParameterData by the parser.
    """
    # suppose at the start that the job is successful
    successful = True
    warnings = []

    if "JOB DONE" not in filecontent:
        successful = False
        msg = "Computation did not finish properly"
        warnings.append(msg)

    if '%%%%%%%%%%%%%%' in filecontent:
        successful = False
        # Pass count=0 to start from the top of the file (anyway, it's a short file)
        # pass an empty warnings list because we don't have existing warnings
        # (this is used to avoid duplication of errors)
        messages = parse_QE_errors(filecontent.split('\n'), count=0,
                                   warnings=[])

        # if it found something, add to log
        if len(messages) > 0:
            warnings.extend(messages)

    out_dict = {
        'warnings': warnings,
    }

    return successful, out_dict


