# -*- coding: utf-8 -*-
"""`Parser` implementation for the `Pw2gwCalculation` calculation job class."""
import io

from aiida.orm import ArrayData, Dict
import numpy as np

from aiida_quantumespresso.utils.mapping import get_logging_container

from .base import BaseParser


class Pw2gwParser(BaseParser):
    """``Parser`` implementation for the ``Pw2gwCalculation`` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed ``Pw2gwCalculation`` into output nodes.

        Two nodes that are expected are the default 'retrieved' `FolderData` node which will store the retrieved files
        permanently in the repository. The second required node is a filepath under the key
        ``retrieved_temporary_files`` which should contain the temporary retrieved files.
        """
        logs = get_logging_container()

        _, parsed_data, logs = self.parse_stdout_from_retrieved(logs)

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        self.out('output_parameters', Dict(dict=parsed_data))

        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE'in logs.error:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE, logs)

        self.exit_code_eps = None
        eps = self.parse_eps_files()

        if self.exit_code_eps:
            return self.exit(self.exit_code_eps, logs)

        self.out('eps', eps)

        return self.exit(logs=logs)

    def parse_eps_files(self):
        """Parse the ``eps*.dat`` files produced by ``pw2gw.x``."""
        retrieved = self.retrieved
        retrieved_names = retrieved.base.repository.list_object_names()

        files = self.node.process_class._internal_retrieve_list
        if any(_ not in retrieved_names for _ in files):
            self.exit_code_eps = self.exit_codes.ERROR_OUTPUT_FILES
            return

        energy = None
        eps = ArrayData()
        for name in self.node.process_class._internal_retrieve_list:
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
