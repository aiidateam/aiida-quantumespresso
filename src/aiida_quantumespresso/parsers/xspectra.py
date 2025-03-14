# -*- coding: utf-8 -*-
import re
from typing import Tuple
import warnings

from aiida.common import AttributeDict
from aiida.orm import Dict, XyData
import numpy as np

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.base import BaseParser
from aiida_quantumespresso.utils.mapping import get_logging_container

warnings.warn(
    'This module is deprecated and will be removed soon as part of migrating XAS and XPS workflows to a new repository.'
    '\nThe new repository can be found at: https://github.com/aiidaplugins/aiida-qe-xspec.', FutureWarning
)

class XspectraParser(BaseParser):
    """Parser for the ``XSpectraCalculation`` calcjob plugin."""

    class_error_map = {
        'Wrong xiabs!!!': 'ERROR_OUTPUT_ABSORBING_SPECIES_WRONG',
        'xiabs < 1 or xiabs > ntyp': 'ERROR_OUTPUT_ABSORBING_SPECIES_ZERO',
        'Calculation not finished': 'ERROR_OUT_OF_WALLTIME',
    }
    success_string='END JOB'

    def parse(self, **kwargs):
        """Parse the contents of the output files stored in the `retrieved` output node."""
        logs = get_logging_container()

        stdout, parsed_data, logs = self.parse_stdout_from_retrieved(logs)

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        parsed_xspectra, logs = self.parse_stdout(stdout, logs)
        parsed_data.update(parsed_xspectra)

        # Parse some additional info which the stdout does not reliably report
        parameters = self.node.inputs.parameters.base.attributes.get('INPUT_XSPECTRA', {})

        xepsilon_defaults = {
            '1': 0,
            '2': 0,
            '3': 1,
        }
        raw_xepsilon = [parameters.get(f'xepsilon({n})', xepsilon_defaults[n]) for n in ['1', '2', '3']]
        parsed_data['xepsilon'] = [float(n) for n in raw_xepsilon]
        parsed_data['xcoordcrys'] = parameters.get('xcoordcrys', True)
        parsed_data['xonly_plot'] = parameters.get('xonly_plot', False)

        self.out('output_parameters', Dict(parsed_data))

        for exit_code in list(self.class_error_map.values()) + ['ERROR_OUTPUT_STDOUT_INCOMPLETE']:
            if exit_code in logs.error:
                return self.exit(self.exit_codes.get(exit_code), logs)

        # Check that the spectra data file exists and is readable
        try:
            with self.retrieved.base.repository.open(self.node.process_class._Spectrum_FILENAME, 'r') as fil:
                xspectra_file = fil.readlines()
        except OSError:
            return self.exit(self.exit_codes.ERROR_READING_SPECTRUM_FILE, logs)

        # Check that the data in the spectra file can be read by NumPy
        try:
            _ = np.genfromtxt(xspectra_file)
        except ValueError:
            return self.exit(self.exit_codes.ERROR_READING_SPECTRUM_FILE_DATA, logs)

        array_names = [[], []]
        array_units = [[], []]

        array_names[0] = ['energy', 'sigma'] # for non-spin-polarised calculations
        array_units[0] = ['eV', 'n/a']

        array_names[1] = ['energy', 'sigma_tot', 'sigma_up', 'sigma_down'] # for spin-polarised calculations
        array_units[1] = ['eV', 'n/a', 'n/a', 'n/a']

        array_data, spin = self.parse_raw_xspectra(xspectra_file, array_names, array_units)

        energy_units = 'eV'
        xy_data = XyData()
        xy_data.set_x(array_data['energy'], 'energy', energy_units)

        y_arrays = []
        y_names = []
        y_units = []
        if spin:
            y_arrays += [array_data['sigma_tot']]
            y_names += ['sigma_tot']
            y_arrays += [array_data['sigma_up']]
            y_arrays += [array_data['sigma_down']]
            y_names += ['sigma_up']
            y_names += ['sigma_down']
            y_units += ['n/a'] * 3
        else:
            y_arrays += [array_data['sigma']]
            y_names += ['sigma']
            y_units += ['n/a']

        xy_data.set_y(y_arrays, y_names, y_units)

        self.out('spectra', xy_data)

        return self.exit(logs=logs)

    @staticmethod
    def parse_stdout(stdout: str, logs: AttributeDict) -> Tuple[dict, AttributeDict]:
        """Parse the ``stdout`` of XSpectra for the core level energy and energy zero of the spectrum."""
        from .parse_raw.base import convert_qe_time_to_sec

        parsed_data = {}

        # Parse the necessary information for data plotting: core level energy of the
        # absorbing atom and the energy zero of the spectrum (typically the Fermi level)
        for line in stdout.split('\n'):
            if 'From SCF save directory' in line:
                if '(spin polarized work)' in line:
                    spin = True
                else:
                    spin = False
                parsed_data['lsda'] = spin
            if 'ehomo [eV]' in line:
                if spin:
                    homo_energy = line.split(':')[-2].split('(')[0].strip()
                else:
                    homo_energy = line.split(':')[-1].split('(')[0].strip()
                homo_energy_units = line.split('[')[1].split(':')[0].replace(']', '')
                parsed_data['highest_occupied_level'] = homo_energy
                parsed_data['highest_occupied_level_units'] = homo_energy_units
            if 'elumo [eV]' in line:
                if spin:
                    lumo_energy = line.split(':')[-2].split('(')[0].strip()
                else:
                    lumo_energy = line.split(':')[-1].split('(')[0].strip()
                lumo_energy_units = line.split('[')[1].split(':')[0].replace(']', '')
                parsed_data['lowest_unoccupied_level'] = lumo_energy
                parsed_data['lowest_unoccupied_level_units'] = lumo_energy_units
                parsed_data['lumo_found'] = True
            elif 'No LUMO value' in line:
                parsed_data['lumo_found'] = False
            if 'ef    [eV]' in line:
                ef_energy = line.split(':')[-1].split('(')[0].strip()
                ef_energy_units = line.split('[')[1].split(':')[0].replace(']', '')
                parsed_data['fermi_energy'] = ef_energy
                parsed_data['fermi_energy_units'] = ef_energy_units

            # parse per-process dynamical RAM estimates
            if 'Estimated max dynamical RAM per process' in line:
                value = line.split('>')[-1]
                match = re.match(r'\s+([+-]?\d+(\.\d*)?|\.\d+([eE][+-]?\d+)?)\s*(Mb|MB|GB)', value)
                if match:
                    try:
                        parsed_data['estimated_ram_per_process'] = float(match.group(1))
                        parsed_data['estimated_ram_per_process_units'] = match.group(4)
                    except (IndexError, ValueError):
                        pass

            # parse total dynamical RAM estimates
            if 'Estimated total dynamical RAM' in line:
                value = line.split('>')[-1]
                match = re.match(r'\s+([+-]?\d+(\.\d*)?|\.\d+([eE][+-]?\d+)?)\s*(Mb|MB|GB)', value)
                if match:
                    try:
                        parsed_data['estimated_ram_total'] = float(match.group(1))
                        parsed_data['estimated_ram_total_units'] = match.group(4)
                    except (IndexError, ValueError):
                        pass

            if 'Core level energy' in line:
                core_energy_line = line
                parsed_data['core_level_energy'] = core_energy_line.split('[')[1].split(':')[1].strip()
                parsed_data['core_level_energy_units'] = core_energy_line.split('[')[1].split(':')[0].replace(']', '')
            if 'energy-zero' in line:
                energy_zero_line = line
                parsed_data['energy_zero'] = energy_zero_line.split('[')[1].split(':')[1].strip()
                parsed_data['energy_zero_units'] = energy_zero_line.split('[')[1].split(':')[0].replace(']', '')

            # Parse the walltime
            # XSpectra does not appear next to the timing data, so we must find 'xanes' instead.
            if 'xanes' in line and 'WALL' in line:
                try:
                    time = line.split('CPU')[1].split('WALL')[0].strip()
                    parsed_data['wall_time'] = time
                except (ValueError, IndexError):
                    break
                try:
                    parsed_data['wall_time_seconds'] = convert_qe_time_to_sec(time)
                except ValueError:
                    logs.warnings.append('Unable to convert wall time from `stdout` to seconds.')

        return parsed_data, logs

    @staticmethod
    def parse_raw_xspectra(xspectra_file, array_names, array_units):
        """Parse the content of the output spectrum.

        This function takes as input the xspectra_file as a list of filelines along with information on how to give
        labels and units to the parsed data.

        :param xspectra_file: xspectra file lines in the form of a list
        :type xspectra_file: list
        :param array_names: list of all array names.
        :type array_names: list
        :param array_units: list of all array units.
        :type array_units: list

        :return array_data: narray, a dictionary for ArrayData type, which
            contains all parsed xspectra output along with labels and units
        """
        xspectra_header = xspectra_file[:4]
        xspectra_data = np.genfromtxt(xspectra_file)

        if len(xspectra_data) == 0:
            raise QEOutputParsingError('XSpectra file is empty.')
        if np.isnan(xspectra_data).any():
            raise QEOutputParsingError('XSpectra file contains non-numeric elements.')

        if len(xspectra_data[0]) == 2:
            array_names = array_names[0]
            array_units = array_units[0]
            spin = False
        elif len(xspectra_data[0]) == 4:
            array_names = array_names[1]
            array_units = array_units[1]
            spin = True
        else:
            raise QEOutputParsingError('XSpectra data file in unsuitable format for the parser')

        i = 0
        array_data = {}
        array_data['header'] = np.array(xspectra_header)
        while i < len(array_names):
            array_data[array_names[i]] = xspectra_data[:, i]
            array_data[array_names[i] + '_units'] = np.array(array_units[i])
            i += 1

        return array_data, spin
