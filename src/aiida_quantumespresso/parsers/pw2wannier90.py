# -*- coding: utf-8 -*-
from aiida.orm import Dict

from aiida_quantumespresso.utils.mapping import get_logging_container

from .base import BaseParser


class Pw2wannier90Parser(BaseParser):
    """``Parser`` implementation for the ``Pw2wannierCalculation`` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed ``Pw2wannierCalculation`` into output nodes.

        Two nodes that are expected are the default 'retrieved' ``FolderData`` node which will store the retrieved files
        permanently in the repository.
        """
        logs = get_logging_container()

        _, parsed_data, logs = self.parse_stdout_from_retrieved(logs)

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        self.out('output_parameters', Dict(parsed_data))

        for exit_code in ['ERROR_OUTPUT_STDOUT_INCOMPLETE']:
            if exit_code in logs.error:
                return self.exit(self.exit_codes.get(exit_code), logs)

        return self.exit(logs=logs)
