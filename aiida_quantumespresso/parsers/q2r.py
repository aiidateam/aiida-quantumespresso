# -*- coding: utf-8 -*-
from __future__ import absolute_import

from aiida.common import exceptions
from aiida.parsers import Parser
from aiida_quantumespresso.calculations.q2r import Q2rCalculation
from aiida_quantumespresso.data.forceconstants import ForceconstantsData


class Q2rParser(Parser):
    """Parser implementation for the Q2rCalculation."""

    def parse(self, **kwargs):
        """Parse the retrieved files from a `Q2rCalculation`."""
        try:
            output_folder = self.retrieved
        except exceptions.NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        filename_stdout = self.node.get_option('output_filename')
        filename_force_constants = Q2rCalculation._FORCE_CONSTANTS_NAME

        if filename_stdout not in output_folder.list_object_names():
            self.logger.error("The standard output file '{}' was not found but is required".format(filename_stdout))
            return self.exit_codes.ERROR_READING_OUTPUT_FILE

        if filename_force_constants not in output_folder.list_object_names():
            self.logger.error("The force constants file '{}' was not found but is required".format(filename_force_constants))
            return self.exit_codes.ERROR_READING_FORCE_CONSTANTS_FILE

        if 'JOB DONE' not in output_folder.get_object_content(filename_stdout):
            self.logger.error('Computation did not finish properly')
            return self.exit_codes.ERROR_JOB_NOT_DONE

        with output_folder.open(filename_force_constants, 'rb') as handle:
            self.out('forceconstants', ForceconstantsData(file=handle))

        return
