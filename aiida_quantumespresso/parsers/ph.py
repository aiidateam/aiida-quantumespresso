# -*- coding: utf-8 -*-
from __future__ import absolute_import

import os
import re
import traceback

from aiida import orm
from aiida.common import exceptions
from aiida.parsers.parser import Parser
from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.parsers.parse_raw.ph import parse_raw_ph_output as parse_stdout


class PhParser(Parser):
    """`Parser` implementation for the `PhCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files from a `PhCalculation`."""
        try:
            self.retrieved
        except exceptions.NotExistent:
            return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_FOLDER)

        # The stdout is required for parsing
        filename_stdout = self.node.get_attribute('output_filename')
        filename_tensor = PhCalculation._OUTPUT_XML_TENSOR_FILE_NAME

        if filename_stdout not in self.retrieved.list_object_names():
            return self.exit_codes.ERROR_OUTPUT_STDOUT_MISSING

        try:
            stdout = self.retrieved.get_object_content(filename_stdout)
        except (IOError, OSError):
            return self.exit_codes.ERROR_OUTPUT_STDOUT_READ

        try:
            tensor_file = self.retrieved.get_object_content(filename_tensor)
        except (IOError, OSError):
            tensor_file = None

        # Look for dynamical matrices
        dynmat_files = []
        dynmat_folder = PhCalculation._FOLDER_DYNAMICAL_MATRIX
        dynmat_prefix = os.path.split(PhCalculation._OUTPUT_DYNAMICAL_MATRIX_PREFIX)[1]

        natural_sort = lambda string: [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', string)]
        for filename in sorted(self.retrieved.list_object_names(dynmat_folder), key=natural_sort):
            if not filename.startswith(dynmat_prefix) or filename.endswith('.freq'):
                continue

            dynmat_files.append(self.retrieved.get_object_content(os.path.join(dynmat_folder, filename)))

        try:
            parsed_data, logs = parse_stdout(stdout, tensor_file, dynmat_files)
        except Exception:
            self.logger.error(traceback.format_exc())
            return self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION

        self.emit_logs(logs)
        self.out('output_parameters', orm.Dict(dict=parsed_data))

        if 'ERROR_OUT_OF_WALLTIME' in logs['error']:
            return self.exit_codes.ERROR_OUT_OF_WALLTIME

        if 'ERROR_CONVERGENCE_NOT_REACHED' in logs['error']:
            return self.exit_codes.ERROR_CONVERGENCE_NOT_REACHED

    def emit_logs(self, *args):
        """Emit the messages in one or multiple "log dictionaries" through the logger of the parser.

        A log dictionary is expected to have the following structure: each key must correspond to a log level of the
        python logging module, e.g. `error` or `warning` and its values must be a list of string messages. The method
        will loop over all log dictionaries and emit the messages it contains with the log level indicated by the key.

        Example log dictionary structure::

            logs = {
                'warning': ['Could not parse the `etot_threshold` variable from the stdout.'],
                'error': ['Self-consistency was not achieved']
            }

        :param args: log dictionaries
        """
        ignore = []

        for logs in args:
            for level, messages in logs.items():
                for message in messages:

                    if message is None:
                        continue

                    stripped = message.strip()

                    if not stripped or stripped in ignore:
                        continue

                    try:
                        getattr(self.logger, level)(stripped)
                    except AttributeError:
                        pass
