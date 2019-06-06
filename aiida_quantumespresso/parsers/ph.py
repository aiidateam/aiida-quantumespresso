# -*- coding: utf-8 -*-
from __future__ import absolute_import

import os

from aiida import orm
from aiida.common import exceptions
from aiida.parsers.parser import Parser
from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.parsers.raw_parser_ph import parse_raw_ph_output


class PhParser(Parser):
    """`Parser` implementation for the `PhCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files from a `PhCalculation`."""
        try:
            output_folder = self.retrieved
        except exceptions.NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        list_of_files = output_folder.list_object_names()

        # The stdout is required for parsing
        filename_stdout = self.node.get_attribute('output_filename')
        filename_tensor = PhCalculation._OUTPUT_XML_TENSOR_FILE_NAME
        filepath_dynmat = PhCalculation._FOLDER_DYNAMICAL_MATRIX

        if filename_stdout not in list_of_files:
            self.logger.error("The standard output file '{}' was not found but is required".format(filename_stdout))
            return self.exit_codes.ERROR_READING_OUTPUT_FILE

        xml_tensor_file = None
        if filename_tensor in list_of_files:
            xml_tensor_file = output_folder._repository._get_base_folder().get_abs_path(filename_tensor)

        # Look for dynamical matrices
        dynmat_dir = output_folder._repository._get_base_folder().get_abs_path(filepath_dynmat)
        dynamical_matrix_list = [
            os.path.join(dynmat_dir, file)
            for file in os.listdir(dynmat_dir)
            if os.path.isfile(os.path.join(dynmat_dir, file)) and
            file.startswith(os.path.split(PhCalculation._OUTPUT_DYNAMICAL_MATRIX_PREFIX)[1]) and
            not file.endswith('.freq')
        ]

        # Sort according to the number at the end of the dynamical matrix files, when there is one
        try:
            dynamical_matrix_list = sorted(dynamical_matrix_list, key=lambda x: int(x.split('-')[-1]))
        except ValueError:
            dynamical_matrix_list = sorted(dynamical_matrix_list)

        out_file = output_folder._repository._get_base_folder().get_abs_path(filename_stdout)
        results, _ = parse_raw_ph_output(out_file, xml_tensor_file, dynamical_matrix_list)

        self.out('output_parameters', orm.Dict(dict=results))

        return
