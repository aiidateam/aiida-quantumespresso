# -*- coding: utf-8 -*-
from __future__ import absolute_import

import numpy as np
from io import StringIO
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
        try:
            retrieved = self.retrieved
        except exceptions.NotExistent:
            return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_FOLDER)
        else:
            retrieved_names = retrieved.list_object_names()

        files = Pw2gwCalculation._internal_retrieve_list + [Pw2gwCalculation._DEFAULT_OUTPUT_FILE]
        if any(not _ in retrieved_names for _ in files):
            return self.exit(self.exit_codes.ERROR_OUTPUT_FILES)

        for name in Pw2gwCalculation._internal_retrieve_list:
            content = retrieved.get_object_content(name)

            try:
                data =  np.loadtxt(StringIO(content))
            except ValueError:
                return self.exit(self.exit_codes.ERROR_OUTPUT_FILES)
            if not len(data.shape) == 2:
                return self.exit(self.exit_codes.ERROR_OUTPUT_FILES)
            if not data.shape[0] > 0 or not data.shape[1] == 2:
                return self.exit(self.exit_codes.ERROR_OUTPUT_FILES)

    def exit(self, exit_code):
        """Log the exit message of the give exit code with level `ERROR` and return the exit code.

        This is a utility function if one wants to return from the parse method and automically add the exit message
        associated to the exit code as a log message to the node: e.g. `return self.exit(self.exit_codes.LABEL))`

        :param exit_code: an `ExitCode`
        :return: the exit code
        """
        self.logger.error(exit_code.message)
        return exit_code



