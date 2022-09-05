# -*- coding: utf-8 -*-
from aiida.orm import Dict

from aiida_quantumespresso.parsers.parse_raw.base import parse_output_base

from .base import Parser


class Pw2wannier90Parser(Parser):
    """`Parser` implementation for the `Pw2wannierCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed `Pw2wannierCalculation` into output nodes.

        Two nodes that are expected are the default 'retrieved' `FolderData` node which will store the retrieved files
        permanently in the repository.
        """
        try:
            filename_stdout = self.node.get_option('output_filename')  # or get_attribute(), but this is clearer
            with self.retrieved.base.repository.open(filename_stdout, 'r') as fil:
                out_file = fil.read()
        except OSError:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_READ)

        parsed_data, logs = parse_output_base(out_file, codename='PW2WANNIER')
        self.emit_logs(logs)

        self.out('output_parameters', Dict(parsed_data))

        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE' in logs.error:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE)
        elif logs.error:
            return self.exit(self.exit_codes.ERROR_GENERIC_QE_ERROR)
