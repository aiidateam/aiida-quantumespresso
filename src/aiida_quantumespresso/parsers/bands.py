# -*- coding: utf-8 -*-
from aiida.orm import Dict

from aiida_quantumespresso.utils.mapping import get_logging_container

from .base import BaseParser


class BandsParser(BaseParser):
    """``Parser`` implementation for the ``BandsCalculation`` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a ``BandsCalculation`` into output nodes."""
        logs = get_logging_container()

        _, parsed_data, logs = self.parse_stdout_from_retrieved(logs)

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        self.out('output_parameters', Dict(parsed_data))

        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE'in logs.error:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE, logs)

        return self.exit(logs=logs)
