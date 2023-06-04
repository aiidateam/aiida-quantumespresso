# -*- coding: utf-8 -*-
"""`Parser` implementation for the `PhCalculation` calculation job class."""
import os
import re

from aiida import orm

from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.parsers.parse_raw.ph import parse_raw_ph_output
from aiida_quantumespresso.utils.mapping import get_logging_container

from .base import BaseParser


class PhParser(BaseParser):
    """``Parser`` implementation for the ``PhCalculation`` calculation job class."""

    class_error_map = {
        'No convergence has been achieved': 'ERROR_CONVERGENCE_NOT_REACHED',
        'problems computing cholesky': 'ERROR_COMPUTING_CHOLESKY',
    }

    def parse(self, **kwargs):
        """Parse the retrieved files from a ``PhCalculation`` into output nodes."""
        logs = get_logging_container()

        stdout, parsed_data, logs = self.parse_stdout_from_retrieved(logs)

        # If the scheduler detected OOW, simply keep that exit code by not returning anything more specific.
        if self.node.exit_status == PhCalculation.exit_codes.ERROR_SCHEDULER_OUT_OF_WALLTIME:
            return

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        filename_tensor = self.node.process_class._OUTPUT_XML_TENSOR_FILE_NAME

        try:
            with self.retrieved.base.repository.open(filename_tensor, 'r') as handle:
                tensor_file = handle.read()
        except OSError:
            tensor_file = None

        # Look for dynamical matrices
        dynmat_files = []
        dynmat_folder = self.node.process_class._FOLDER_DYNAMICAL_MATRIX
        dynmat_prefix = os.path.split(self.node.process_class._OUTPUT_DYNAMICAL_MATRIX_PREFIX)[1]

        natural_sort = lambda string: [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', string)]
        for filename in sorted(self.retrieved.base.repository.list_object_names(dynmat_folder), key=natural_sort):
            if not filename.startswith(dynmat_prefix) or filename.endswith('.freq'):
                continue

            dynmat_files.append(self.retrieved.base.repository.get_object_content(os.path.join(dynmat_folder, filename)))

        parsed_ph_data, logs = parse_raw_ph_output(stdout, logs, tensor_file, dynmat_files)
        parsed_data.update(parsed_ph_data)

        self.out('output_parameters', orm.Dict(parsed_data))

        for exit_code in list(self.get_error_map().values()) + ['ERROR_OUTPUT_STDOUT_INCOMPLETE']:
            if exit_code in logs.error:
                return self.exit(self.exit_codes.get(exit_code), logs)

        return self.exit(logs=logs)
