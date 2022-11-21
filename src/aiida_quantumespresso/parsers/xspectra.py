# -*- coding: utf-8 -*-
import re

from aiida.orm import Dict, XyData
import numpy as np

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.base import Parser


class XspectraParser(Parser):
    """ Parser for the XSpectraCalculation calcjob plugin """

    def parse(self, **kwargs):
        """Parse the contents of the output files stored in the `retrieved` output node."""
        from aiida.plugins import DataFactory

        retrieved = self.retrieved
        try:
            filename_stdout = self.node.get_option('output_filename')
            with retrieved.base.repository.open(filename_stdout, 'r') as fil:
                out_file = fil.readlines()
        except OSError:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_READ)

        # Check the stdout for obvious errors
        job_done = False
        for line in out_file:
            if 'Wrong xiabs!!!' in line:
                return self.exit(self.exit_codes.ERROR_OUTPUT_ABSORBING_SPECIES_WRONG)
            if 'xiabs < 1 or xiabs > ntyp' in line:
                return self.exit(self.exit_codes.ERROR_OUTPUT_ABSORBING_SPECIES_ZERO)
            if 'Calculation not finished' in line:
                return self.exit(self.exit_codes.ERROR_OUT_OF_WALLTIME)
            if 'END JOB' in line:
                job_done = True
                break
        if not job_done:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE)

        # Check that the spectra data file exists and is readable
        try:
            with retrieved.base.repository.open(self.node.process_class._Spectrum_FILENAME, 'r') as fil:
                xspectra_file = fil.readlines()
        except OSError:
            return self.exit(self.exit_codes.ERROR_READING_SPECTRUM_FILE)

        # Check that the data in the spectra file can be read by NumPy
        try:
            xspectra_data = np.genfromtxt(xspectra_file)
        except ValueError:
            return self.exit(self.exit_codes.ERROR_READING_SPECTRUM_FILE_DATA)

        # end of initial checks

        array_names = [[], []]
        array_units = [[], []]

        array_names[0] = ['energy', 'sigma'] # for non-spin-polarised calculations
        array_units[0] = ['eV', 'n/a']

        array_names[1] = ['energy', 'sigma_tot', 'sigma_up', 'sigma_down'] # for spin-polarised calculations
        array_units[1] = ['eV', 'n/a', 'n/a', 'n/a']

        array_data, spin = parse_raw_xspectra(xspectra_file, array_names, array_units)

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

        parsed_data, logs = parse_stdout_xspectra(filecontent=out_file, codename='XSpectra')

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
        self.emit_logs(logs)

        self.out('spectra', xy_data)
        self.out('output_parameters', Dict(dict=parsed_data))

def parse_raw_xspectra(xspectra_file, array_names, array_units):
    """Parse the content of the output spectrum.

    This function takes as input the xspectra_file as a list of filelines along with information on how to give labels
    and units to the parsed data.

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

def parse_stdout_xspectra(filecontent, codename=None, message_map=None):
    """Parses the output file of an XSpectra calculation, checking for
    basic content like END JOB, errors with %%%%, and the core level energy
    and the energy zero of the spectrum.

    :param filecontent: a string with the output file content
    :param codename: the string printed both in the header and near the
                    walltime. If passed, a few more things are parsed (e.g.
                    code version, walltime, ...)
    :returns: tuple of two dictionaries, with the parsed data and log
              messages, respectively
    """
    from aiida_quantumespresso.utils.mapping import get_logging_container

    from .parse_raw.base import convert_qe_time_to_sec

    keys = ['error', 'warning']

    if message_map is not None and (not isinstance(message_map, dict) or any(key not in message_map for key in keys)):
        raise RuntimeError(f'invalid format `message_map`: should be dictionary with two keys {keys}')

    logs = get_logging_container()
    parsed_data = {}

    lines = filecontent if isinstance(filecontent, list) else filecontent.split('\n')

    # Parse the necessary information for data plotting: core level energy of the
    # absorbing atom and the energy zero of the spectrum (typically the Fermi level)
    for line in lines:
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

    if codename is not None:

        codestring = f'Program {codename}'

        for line_number, line in enumerate(lines):

            if codestring in line and 'starts on' in line:
                parsed_data['code_version'] = line.split(codestring)[1].split('starts on')[0].strip()

            # Parse the walltime
            # XSpectra does not appear next to the timing data, so we must find 'xanes' instead.
            if 'xanes' in line and 'WALL' in line:
                    try:
                        time = line.split('CPU')[1].split('WALL')[0].strip()
                        parsed_data['wall_time'] = time
                    except (ValueError, IndexError):
                        logs.warnings.append('ERROR_PARSING_WALLTIME')
                    try:
                        parsed_data['wall_time_seconds'] = convert_qe_time_to_sec(time)
                    except ValueError:
                        logs.warnings.append('ERROR_CONVERTING_WALLTIME_TO_SECONDS')

            # Parse an error message with optional mapping of the message
            if '%%%%%%%%%%%%%%' in line:
                parse_output_error(lines, line_number, logs, message_map)

    return parsed_data, logs
