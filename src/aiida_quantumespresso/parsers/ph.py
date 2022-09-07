# -*- coding: utf-8 -*-
"""`Parser` implementation for the `PhCalculation` calculation job class."""
import os
import re
import traceback

from aiida import orm

from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.parsers.parse_raw.ph import parse_raw_ph_output as parse_stdout

from .base import Parser


class PhParser(Parser):
    """`Parser` implementation for the `PhCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files from a `PhCalculation`."""
        retrieved = self.retrieved

        # The stdout is required for parsing
        filename_stdout = self.node.base.attributes.get('output_filename')
        filename_tensor = PhCalculation._OUTPUT_XML_TENSOR_FILE_NAME

        if filename_stdout not in retrieved.base.repository.list_object_names():
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_MISSING)

        try:
            stdout = retrieved.base.repository.get_object_content(filename_stdout)
        except (IOError, OSError):
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_READ)

        try:
            tensor_file = retrieved.base.repository.get_object_content(filename_tensor)
        except (IOError, OSError):
            tensor_file = None

        # Look for dynamical matrices
        dynmat_files = []
        dynmat_folder = PhCalculation._FOLDER_DYNAMICAL_MATRIX  # pylint: disable=protected-access
        dynmat_prefix = os.path.split(PhCalculation._OUTPUT_DYNAMICAL_MATRIX_PREFIX)[1]  # pylint: disable=protected-access

        natural_sort = lambda string: [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', string)]
        for filename in sorted(retrieved.base.repository.list_object_names(dynmat_folder), key=natural_sort):
            if not filename.startswith(dynmat_prefix) or filename.endswith('.freq'):
                continue

            dynmat_files.append(retrieved.base.repository.get_object_content(os.path.join(dynmat_folder, filename)))

        try:
            parsed_data, logs = parse_stdout(stdout, tensor_file, dynmat_files)
        except Exception as exc:
            self.logger.error(traceback.format_exc())
            return self.exit(self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION.format(exception=exc))

        self.emit_logs(logs)
        self.out('output_parameters', orm.Dict(parsed_data))

        # If the scheduler detected OOW, simply keep that exit code by not returning anything more specific.
        if self.node.exit_status == PhCalculation.exit_codes.ERROR_SCHEDULER_OUT_OF_WALLTIME:
            return

        if 'ERROR_OUT_OF_WALLTIME' in logs['error']:
            return self.exit_codes.ERROR_OUT_OF_WALLTIME

        if 'ERROR_CONVERGENCE_NOT_REACHED' in logs['error']:
            return self.exit_codes.ERROR_CONVERGENCE_NOT_REACHED

        if 'ERROR_COMPUTING_CHOLESKY' in logs['error']:
            return self.exit_codes.ERROR_COMPUTING_CHOLESKY

        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE' in logs['error']:
            return self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE
