# -*- coding: utf-8 -*-
from __future__ import absolute_import

from aiida.common import NotExistent
from aiida.orm import Dict
from aiida.parsers import Parser
from aiida_quantumespresso.parsers.parse_raw.base import parse_output_base, emit_logs


class Pw2wannier90Parser(Parser):
    """This class is the implementation of the Parser class for pw2wannier90.x."""

    def parse(self, **kwargs):
        """Parses the datafolder, stores results.

        In this case we only parse the aiida.out file, and retrieve any files given in the internal and additional
        retrieve lists.
        """
        # Check that the retrieved folder is there
        try:
            out_folder = self.retrieved
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # Read standard out
        try:
            filename_stdout = self.node.get_option('output_filename')  # or get_attribute(), but this is clearer
            with out_folder.open(filename_stdout, 'r') as fil:
                out_file = fil.read()
        except OSError:
            return self.exit_codes.ERROR_OUTPUT_STDOUT_READ

        # check that the file has finished (i.e. JOB DONE is inside the file)
        out_dict, logs = parse_output_base(out_file, codename='PW2WANNIER')
        emit_logs(self.logger, logs)

        # Output a Dict with whatever has been parsed
        self.out('output_parameters', Dict(dict=out_dict))

        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE' in logs.error:
            return self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE
        elif logs.error:
            return self.exit_codes.ERROR_GENERIC_QE_ERROR
