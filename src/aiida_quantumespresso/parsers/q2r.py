# -*- coding: utf-8 -*-
from aiida.orm import Dict

from aiida_quantumespresso.data.force_constants import ForceConstantsData

from .base import BaseParser


class Q2rParser(BaseParser):
    """``Parser`` implementation for the ``Q2rCalculation`` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a ``Q2rCalculation`` into output nodes."""
        parsed_stdout, logs_stdout = self._parse_stdout_from_retrieved()
        self._emit_logs(logs_stdout)

        self.out('output_parameters', Dict(parsed_stdout))

        for exit_code in ['ERROR_OUTPUT_STDOUT_MISSING', 'ERROR_OUTPUT_STDOUT_READ', 'ERROR_OUTPUT_STDOUT_INCOMPLETE']:
            if exit_code in logs_stdout.error:
                return self._exit(self.exit_codes.get(exit_code))

        filename_force_constants = self.node.process_class._FORCE_CONSTANTS_NAME

        if filename_force_constants not in self.retrieved.base.repository.list_object_names():
            return self._exit(self.exit_codes.ERROR_READING_FORCE_CONSTANTS_FILE)

        with self.retrieved.base.repository.open(filename_force_constants, 'rb') as handle:
            self.out('force_constants', ForceConstantsData(file=handle))

        return
