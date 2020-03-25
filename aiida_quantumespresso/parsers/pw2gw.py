# -*- coding: utf-8 -*-
from __future__ import absolute_import

import numpy as np
from io import StringIO
from aiida import orm
from aiida.common import exceptions
from aiida.parsers import Parser

from aiida_quantumespresso.calculations.pw2gw import Pw2gwCalculation

class Pw2gwParser(Parser):
    """`Parser` implementation for the `Pw2gwCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed `Pw2gwCalculation` into output nodes.

        Two nodes that are expected are the default 'retrieved' `FolderData` node which will store the retrieved files
        permanently in the repository. The second required node is a filepath under the key `retrieved_temporary_files`
        which should contain the temporary retrieved files.
        """
        self.exit_code_stdout = None
        self.exit_code_eps    = None

        try:
            self.retrieved
        except exceptions.NotExistent:
            return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_FOLDER)

        # Parse the pw2gw stout file
        data, logs_stdout = self.parse_stdout()

        self.emit_logs(logs_stdout)

        if self.exit_code_stdout:
            return self.exit(self.exit_code_stdout)

        self.out('output_parameters', orm.Dict(dict=data))

        # Parse the pw2g outputfiles
        eps = self.parse_eps_files()

        if self.exit_code_eps:
            return self.exit(self.exit_code_eps)

        self.out('eps', eps)


    def parse_eps_files(self):
        """Parse the eps*.dat files produced by pw2gw.x and store them in the `eps` node."""
        retrieved = self.retrieved
        retrieved_names = retrieved.list_object_names()

        files = Pw2gwCalculation._internal_retrieve_list
        if any(not _ in retrieved_names for _ in files):
            self.exit_code_eps = self.exit_codes.ERROR_OUTPUT_FILES
            return

        energy = None
        eps = orm.ArrayData()
        for name in Pw2gwCalculation._internal_retrieve_list:
            content = retrieved.get_object_content(name)
            base = name.split('.')[0]

            try:
                data =  np.loadtxt(StringIO(content))
            except ValueError:
                self.exit_code_eps = self.exit_codes.ERROR_OUTPUT_FILES
                return
            if len(data.shape) != 2 or data.shape[0] == 0 or data.shape[1] != 2:
                self.exit_code_eps = self.exit_codes.ERROR_OUTPUT_FILES_INVALID_FORMAT
                return

            X, Y = data.T
            if energy is None:
                energy = X
                eps.set_array('energy', X)
            elif not np.allclose(X, energy):
                self.exit_code_eps = self.exit_codes.ERROR_OUTPUT_FILES_ENERGY_MISMATCH
                return

            eps.set_array(base, Y)

        return eps

    def parse_stdout(self):
        """Parse the stdout file of pw2gw to build the `output_parameters` node."""
        from aiida_quantumespresso.utils.mapping import get_logging_container
        from aiida_quantumespresso.parsers.parse_raw.pw2gw import parse_stdout

        logs = get_logging_container()
        parsed_data = {}

        filename_stdout = self.node.get_attribute('output_filename')
        # print(" < ", filename_stdout, ' > ')
        # return parsed_data, logs

        if filename_stdout not in self.retrieved.list_object_names():
            self.exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_MISSING
            return parsed_data, logs

        try:
            stdout = self.retrieved.get_object_content(filename_stdout)
        except IOError:
            self.exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_READ
            return parsed_data, logs

        try:
            parsed_data, logs = parse_stdout(stdout)
        except Exception:
            import traceback
            traceback.print_exc()
            self.exit_code_stdout = self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION

        # If the stdout was incomplete, most likely the job was interrupted before it could cleanly finish, so the
        # output files are most likely corrupt and cannot be restarted from
        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE' in logs['error']:
            self.exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE

        return parsed_data, logs

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
        ignore = [
            'Error while parsing ethr.',
            'DEPRECATED: symmetry with ibrav=0, use correct ibrav instead'
        ]

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

    def exit(self, exit_code):
        """Log the exit message of the give exit code with level `ERROR` and return the exit code.

        This is a utility function if one wants to return from the parse method and automically add the exit message
        associated to the exit code as a log message to the node: e.g. `return self.exit(self.exit_codes.LABEL))`

        :param exit_code: an `ExitCode`
        :return: the exit code
        """
        self.logger.error(exit_code.message)
        return exit_code



