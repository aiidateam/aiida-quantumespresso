# -*- coding: utf-8 -*-
"""`Parser` implementation for the `Pw2gwCalculation` calculation job class."""
import io

from aiida import orm
import numpy as np

from aiida_quantumespresso.calculations.pw2gw import Pw2gwCalculation

from .base import Parser


class Pw2gwParser(Parser):
    """`Parser` implementation for the `Pw2gwCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed `Pw2gwCalculation` into output nodes.

        Two nodes that are expected are the default 'retrieved' `FolderData` node which will store the retrieved files
        permanently in the repository. The second required node is a filepath under the key `retrieved_temporary_files`
        which should contain the temporary retrieved files.
        """
        self.exit_code_stdout = None
        self.exit_code_eps = None

        # Parse the pw2gw stout file
        data, logs_stdout = self.parse_stdout()

        self.emit_logs(logs_stdout)

        if self.exit_code_stdout:
            return self.exit(self.exit_code_stdout)

        self.out('output_parameters', orm.Dict(data))

        # Parse the pw2g outputfiles
        eps = self.parse_eps_files()

        if self.exit_code_eps:
            return self.exit(self.exit_code_eps)

        self.out('eps', eps)

    def parse_eps_files(self):
        """Parse the eps*.dat files produced by pw2gw.x and store them in the `eps` node."""
        retrieved = self.retrieved
        retrieved_names = retrieved.base.repository.list_object_names()

        files = Pw2gwCalculation._internal_retrieve_list
        if any(_ not in retrieved_names for _ in files):
            self.exit_code_eps = self.exit_codes.ERROR_OUTPUT_FILES
            return

        energy = None
        eps = orm.ArrayData()
        for name in Pw2gwCalculation._internal_retrieve_list:
            content = retrieved.base.repository.get_object_content(name)
            base = name.split('.')[0]

            try:
                data = np.loadtxt(io.StringIO(content))
            except ValueError:
                self.exit_code_eps = self.exit_codes.ERROR_OUTPUT_FILES
                return
            if len(data.shape) != 2 or data.shape[0] == 0 or data.shape[1] != 2:
                self.exit_code_eps = self.exit_codes.ERROR_OUTPUT_FILES_INVALID_FORMAT
                return

            x, y = data.T
            if energy is None:
                energy = x
                eps.set_array('energy', x)
            elif not np.allclose(x, energy):
                self.exit_code_eps = self.exit_codes.ERROR_OUTPUT_FILES_ENERGY_MISMATCH
                return

            eps.set_array(base, y)

        return eps

    def parse_stdout(self):
        """Parse the stdout file of pw2gw to build the `output_parameters` node."""
        from aiida_quantumespresso.parsers.parse_raw.pw2gw import parse_stdout
        from aiida_quantumespresso.utils.mapping import get_logging_container

        logs = get_logging_container()
        parsed_data = {}

        filename_stdout = self.node.base.attributes.get('output_filename')

        if filename_stdout not in self.retrieved.base.repository.list_object_names():
            self.exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_MISSING
            return parsed_data, logs

        try:
            stdout = self.retrieved.base.repository.get_object_content(filename_stdout)
        except IOError:
            self.exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_READ
            return parsed_data, logs

        try:
            parsed_data, logs = parse_stdout(stdout)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            self.exit_code_stdout = self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION.format(exception=exc)

        # If the stdout was incomplete, most likely the job was interrupted before it could cleanly finish, so the
        # output files are most likely corrupt and cannot be restarted from
        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE' in logs['error']:
            self.exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE

        return parsed_data, logs
