# -*- coding: utf-8 -*-
from aiida.orm import Dict

from .base import BaseParser


class Pw2wannier90Parser(BaseParser):
    """``Parser`` implementation for the ``Pw2wannierCalculation`` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed ``Pw2wannierCalculation`` into output nodes.

        Two nodes that are expected are the default 'retrieved' ``FolderData`` node which will store the retrieved files
        permanently in the repository.
        """
        parsed_stdout, logs_stdout = self._parse_stdout_from_retrieved()
        self._emit_logs(logs_stdout)

        self.out('output_parameters', Dict(parsed_stdout))

        for exit_code in ['ERROR_OUTPUT_STDOUT_MISSING', 'ERROR_OUTPUT_STDOUT_READ', 'ERROR_OUTPUT_STDOUT_INCOMPLETE']:
            if exit_code in logs_stdout.error:
                return self._exit(self.exit_codes.get(exit_code))
