# -*- coding: utf-8 -*-
from __future__ import absolute_import
from aiida.parsers.parser import Parser
from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.calculations.q2r import Q2rCalculation


class Q2rParser(Parser):
    """Parser implementation for the Q2rCalculation."""

    def __init__(self, calculation):
        if not isinstance(calculation, Q2rCalculation):
            raise QEOutputParsingError('Input calc must be a Q2rCalculation')

        self._calc = calculation

        super(Q2rParser, self).__init__(calculation)

    def parse(self, **kwargs):
        """Parse the output data stored in the retrieved folder."""
        successful = True
        nodes = []

        # Check that the retrieved folder is there
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
        except KeyError:
            self.logger.error('No retrieved folder found')
            return False, ()

        # At least the stdout should exist
        if self._calc._OUTPUT_FILE_NAME not in out_folder.get_folder_list():
            self.logger.error('Standard output not found')
            return False, ()

        # Check that the file has finished (i.e. JOB DONE is inside the file)
        filpath = out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME)
        with open(filpath, 'r') as handle:
            lines = handle.read()

        if 'JOB DONE' not in lines:
            successful = False
            self.logger.error('Computation did not finish properly')

        # Check that there is one and only one real space force constant matrix present
        try:
            link_name = self._calc.get_linkname_force_matrix()
            self._calc.get_outgoing(link_label_filter=link_name).one()
        except ValueError:
            successful = False
            self.logger.error('Found no or more than one output nodes with the link label ()'.format(link_name))

        return successful, nodes
