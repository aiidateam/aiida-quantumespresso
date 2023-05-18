# -*- coding: utf-8 -*-
from aiida.orm import Dict, XyData
import numpy as np

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.utils.mapping import get_logging_container

from .base import BaseParser


class DosParser(BaseParser):
    """``Parser`` implementation for the ``DosCalculation`` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a ``DosCalculation`` into output nodes."""
        logs = get_logging_container()

        _, parsed_stdout, logs = self.parse_stdout_from_retrieved(logs)

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        self.out('output_parameters', Dict(parsed_stdout))

        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE'in logs.error:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE, logs)

        # Parse the DOS
        try:
            with self.retrieved.base.repository.open(self.node.process_class._DOS_FILENAME, 'r') as handle:
                dos_file = handle.readlines()
        except OSError:
            return self.exit(self.exit_codes.ERROR_READING_DOS_FILE, logs)

        array_names = [[], []]
        array_units = [[], []]
        array_names[0] = ['dos_energy', 'dos', 'integrated_dos']  # When spin is not displayed
        array_names[1] = ['dos_energy', 'dos_spin_up', 'dos_spin_down', 'integrated_dos']  # When spin is displayed
        array_units[0] = ['eV', 'states/eV', 'states']  # When spin is not displayed
        array_units[1] = ['eV', 'states/eV', 'states/eV', 'states']  # When spin is displayed

        # grabs parsed data from aiida.dos
        # TODO: should I catch any QEOutputParsingError from parse_raw_dos,
        #       log an error and return an exit code?
        array_data, spin = parse_raw_dos(dos_file, array_names, array_units)

        energy_units = 'eV'
        dos_units = 'states/eV'
        int_dos_units = 'states'
        xy_data = XyData()
        xy_data.set_x(array_data['dos_energy'], 'dos_energy', energy_units)
        y_arrays = []
        y_names = []
        y_units = []
        y_arrays += [array_data['integrated_dos']]
        y_names += ['integrated_dos']
        y_units += ['states']
        if spin:
            y_arrays += [array_data['dos_spin_up']]
            y_arrays += [array_data['dos_spin_down']]
            y_names += ['dos_spin_up']
            y_names += ['dos_spin_down']
            y_units += ['states/eV'] * 2
        else:
            y_arrays += [array_data['dos']]
            y_names += ['dos']
            y_units += ['states/eV']
        xy_data.set_y(y_arrays, y_names, y_units)

        self.out('output_dos', xy_data)

        return self.exit(logs=logs)


def parse_raw_dos(dos_file, array_names, array_units):
    """This function takes as input the dos_file as a list of filelines along with information on how to give labels and
    units to the parsed data.

    :param dos_file: dos file lines in the form of a list
    :type dos_file: list
    :param array_names: list of all array names, note that array_names[0]
                        is for the case with non spin-polarized calculations
                        and array_names[1] is for the case with spin-polarized
                        calculation
    :type array_names: list
    :param array_units: list of all array units, note that array_units[0] is
                        for the case with non spin-polarized calculations and
                        array_units[1] is for the case with spin-polarized
                        calculation
    :type array_units: list

    :return array_data: narray, a dictionary for ArrayData type, which contains
                        all parsed dos output along with labels and units
    :return spin: boolean, indicates whether the parsed results are spin
                  polarized
    """

    dos_header = dos_file[0]
    try:
        dos_data = np.genfromtxt(dos_file)
    except ValueError:
        raise QEOutputParsingError('dosfile could not be loaded using genfromtxt')
    if len(dos_data) == 0:
        raise QEOutputParsingError('Dos file is empty.')
    if np.isnan(dos_data).any():
        raise QEOutputParsingError('Dos file contains non-numeric elements.')

    # Checks the number of columns, essentially to see whether spin was used
    if len(dos_data[0]) == 3:
        # spin is not used
        array_names = array_names[0]
        array_units = array_units[0]
        spin = False
    elif len(dos_data[0]) == 4:
        # spin is used
        array_names = array_names[1]
        array_units = array_units[1]
        spin = True
    else:
        raise QEOutputParsingError('Dos file in format that the parser is not designed to handle.')

    i = 0
    array_data = {}
    array_data['header'] = np.array(dos_header)
    while i < len(array_names):
        array_data[array_names[i]] = dos_data[:, i]
        array_data[array_names[i] + '_units'] = np.array(array_units[i])
        i += 1
    return array_data, spin
