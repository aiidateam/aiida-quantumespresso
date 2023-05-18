# -*- coding: utf-8 -*-
from aiida.orm import Dict

from aiida_quantumespresso.data.force_constants import ForceConstantsData
from aiida_quantumespresso.utils.mapping import get_logging_container

from .base import BaseParser


class Q2rParser(BaseParser):
    """``Parser`` implementation for the ``Q2rCalculation`` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a ``Q2rCalculation`` into output nodes."""
        logs = get_logging_container()

        _, parsed_data, logs = self.parse_stdout_from_retrieved(logs)

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        self.out('output_parameters', Dict(parsed_data))

        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE'in logs.error:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE, logs)

        filename_force_constants = self.node.process_class._FORCE_CONSTANTS_NAME

        if filename_force_constants not in self.retrieved.base.repository.list_object_names():
            return self.exit(self.exit_codes.ERROR_READING_FORCE_CONSTANTS_FILE, logs)

        with self.retrieved.base.repository.open(filename_force_constants, 'rb') as handle:
            self.out('force_constants', ForceConstantsData(file=handle))

        return self.exit(logs=logs)
