# -*- coding: utf-8 -*-
"""`Parser` implementation for the `PhCalculation` calculation job class."""
import os
import re

from aiida import orm

from aiida_quantumespresso.calculations.ph import PhCalculation

from .base import BaseParser


class PhParser(BaseParser):
    """``Parser`` implementation for the ``PhCalculation`` calculation job class."""

    class_error_map = {
        'No convergence has been achieved': 'ERROR_CONVERGENCE_NOT_REACHED',
        'problems computing cholesky': 'ERROR_COMPUTING_CHOLESKY',
    }

    def parse(self, **kwargs):
        """Parse the retrieved files from a ``PhCalculation`` into output nodes."""
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

        parsed_stdout, logs_stdout = self._parse_stdout_from_retrieved(
            tensor_file=tensor_file, dynmat_files=dynmat_files
        )
        self._emit_logs(logs_stdout)
        self.out('output_parameters', orm.Dict(parsed_stdout))

        for exit_code in [
            'ERROR_OUTPUT_STDOUT_MISSING', 'ERROR_OUTPUT_STDOUT_READ', 'ERROR_OUT_OF_WALLTIME',
            'ERROR_CONVERGENCE_NOT_REACHED', 'ERROR_COMPUTING_CHOLESKY', 'ERROR_OUTPUT_STDOUT_INCOMPLETE'
        ]:
            if exit_code in logs_stdout.error:
                return self._exit(self.exit_codes.get(exit_code))

        # If the scheduler detected OOW, simply keep that exit code by not returning anything more specific.
        if self.node.exit_status == PhCalculation.exit_codes.ERROR_SCHEDULER_OUT_OF_WALLTIME:
            return

    @classmethod
    def parse_stdout(cls, stdout: str, tensor_file: str, dynmat_files: list) -> tuple:
        """Parse the ``stdout`` content of a Quantum ESPRESSO ``ph.x`` calculation.

        :param stdout: the stdout content as a string.
        :returns: tuple of two dictionaries, with the parsed data and log messages, respectively.
        """
        from aiida_quantumespresso.parsers.parse_raw.ph import parse_raw_ph_output as parse_stdout

        parsed_base, logs_base = super().parse_stdout(stdout)

        # if logs_base['error']:
        #     return parsed_base, logs_base

        parsed_data, logs = parse_stdout(stdout, tensor_file, dynmat_files)



        parsed_data.update(parsed_base)
        for log_level, log_items in logs_base.items():
            logs[log_level].extend(log_items)

        return parsed_data, logs
