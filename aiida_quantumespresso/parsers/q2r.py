# -*- coding: utf-8 -*-
from aiida.orm.data.folder import FolderData
from aiida.parsers.parser import Parser
from aiida.common.datastructures import calc_states
from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.calculations.q2r import Q2rCalculation


class Q2rParser(Parser):
    """
    This class is the implementation of the Parser class for Q2r.
    """

    def __init__(self, calculation):
        """
        Initialize the instance of Q2rParser
        """
        # check for valid input
        if not isinstance(calculation, Q2rCalculation):
            raise QEOutputParsingError("Input calc must be a Q2rCalculation")

        self._calc = calculation

        super(Q2rParser, self).__init__(calculation)

    def parse_with_retrieved(self, retrieved):
        """      
        Parses the datafolder, stores results.
        This parser for this simple code does simply store in the DB a node
        representing the file of forces in real space
        """
        from aiida.common.exceptions import InvalidOperation

        # suppose at the start that the job is successful
        successful = True

        # Check that the retrieved folder is there
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
        except KeyError:
            self.logger.error("No retrieved folder found")
            return False, ()

        # check what is inside the folder
        list_of_files = out_folder.get_folder_list()
        # at least the stdout should exist
        if not self._calc._OUTPUT_FILE_NAME in list_of_files:
            successful = False
            self.logger.error("Standard output not found")
            return successful, ()

        # check that the file has finished (i.e. JOB DONE is inside the file)
        filpath = out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME)
        with open(filpath, 'r') as fil:
            lines = fil.read()
        if "JOB DONE" not in lines:
            successful = False
            self.logger.error("Computation did not finish properly")

        # check that the real space force constant matrix is present
        the_outputs = self._calc.get_outputs(also_labels=True)
        the_files = [i[1] for i in the_outputs if i[0] == self._calc.get_linkname_force_matrix()]

        if len(the_files) != 1:
            successful = False
            self.logger.error("There should be only one force constants file; "
                              "found instead {} file(s)".format(len(the_files)))

        new_nodes_list = []

        return successful, new_nodes_list
