# -*- coding: utf-8 -*-
from aiida.common import exceptions

from aiida_quantumespresso.calculations.q2r import Q2rCalculation
from aiida_quantumespresso.data.force_constants import ForceConstantsData
from .base import Parser


class Q2rParser(Parser):
    """Parser implementation for the Q2rCalculation."""

    def parse(self, **kwargs):
        """Parse the retrieved files from a `Q2rCalculation`."""
        try:
            output_folder = self.retrieved
        except exceptions.NotExistent:
            return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_FOLDER)

        filename_stdout = self.node.get_option('output_filename')
        filename_force_constants = Q2rCalculation._FORCE_CONSTANTS_NAME

        if filename_stdout not in output_folder.list_object_names():
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_READ)

        if filename_force_constants not in output_folder.list_object_names():
            return self.exit(self.exit_codes.ERROR_READING_FORCE_CONSTANTS_FILE)

        if 'JOB DONE' not in output_folder.get_object_content(filename_stdout):
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE)

        with output_folder.open(filename_force_constants, 'rb') as handle:
            self.out('force_constants', ForceConstantsData(file=handle))

        return
