# -*- coding: utf-8 -*-
"""`Parser` implementation for the `average.x` code of Quantum ESPRESSO."""

import os

from aiida import orm
from aiida.plugins import CalculationFactory
from ase.units import Bohr, Ry
import numpy as np

from aiida_quantumespresso.calculations.average import AverageCalculation
from aiida_quantumespresso.utils.mapping import get_logging_container

from .base import BaseParser

AverageCalculation = CalculationFactory("average")


class AverageParser(BaseParser):
    """
    Parser for the output of average.x.

    This parser produces two outputs:
    - 'output_parameters': a Dict containing parsed metadata from the standard output.
    - 'output_data': an ArrayData node with a single array named 'average', containing three columns: z (Angstrom), p(z) (eV), and m(z) (eV).

    Assumes the presence of a data file named as specified by AverageCalculation._DEFAULT_OUTPUT_DATA_FILE
    in the retrieved temporary folder, formatted as columns of z (Bohr), p(z) (Ry), and m(z) (Ry).
    """

    def parse(self, **kwargs):
        """Parse the retrieved files into output nodes."""
        logs = get_logging_container()

        _, parsed_data, logs = self.parse_stdout_from_retrieved(logs)

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        self.out("output_parameters", orm.Dict(parsed_data))

        try:
            retrieved_temporary_folder = kwargs["retrieved_temporary_folder"]
        except KeyError:
            return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_TEMPORARY_FOLDER)

        # pylint: disable=protected-access
        try:
            with open(
                os.path.join(
                    retrieved_temporary_folder,
                    AverageCalculation._DEFAULT_OUTPUT_DATA_FILE,
                ),
                "r",
                encoding="utf-8",
            ) as file_handler:
                self.out("output_data", self.parse_data(file_handler))
        except OSError:
            return self.exit(
                self.exit_codes.ERROR_OUTPUT_DATAFILE_READ.format(
                    filename=AverageCalculation._DEFAULT_OUTPUT_DATA_FILE
                ),
                logs,
            )
        except Exception as exception:  # pylint: disable=broad-except
            return self.exit(
                self.exit_codes.ERROR_OUTPUT_DATAFILE_PARSE.format(
                    filename=AverageCalculation._DEFAULT_OUTPUT_DATA_FILE,
                    exception=exception,
                ),
                logs,
            )

    def parse_data(self, file_handler):
        """
        Parse the data file.

        The data file contains the following columns:
        - z Bohr
        - p(z) Rydberg
        - m(z) Rydberg

        Output ArrayData node will contain one array with three columns:
        - z Ang
        - p(z) eV
        - m(z) eV
        """
        data = np.loadtxt(file_handler)
        data[:, 0] *= Bohr
        data[:, 1:] *= Ry

        return orm.ArrayData(data)
