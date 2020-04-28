# -*- coding: utf-8 -*-
"""`Parser` implementation for the `PpCalculation` calculation job class."""
from __future__ import absolute_import

import traceback
import os.path

import numpy as np

from aiida import orm
from aiida.common import exceptions

from aiida_quantumespresso.calculations.pp import PpCalculation
from aiida_quantumespresso.utils.mapping import get_logging_container
from .base import Parser

from six.moves import range
from six.moves import zip


class PpParser(Parser):
    """`Parser` implementation for the `PpCalculation` calculation job class."""

    # Lookup: plot_num --> units
    units_dict = {
        0: 'e/bohr^3',  # Electrons, electronic charge density
        1: 'Ry',  # Total potential
        2: 'Ry',  # Ionic potential
        3: 'states/bohr^3',  # Density of states over an energy range
        4: 'Ry/K.bohr^3',  # Local density of electronic entropy
        5: 'states/bohr^3',  # Simulated STM images from LDOS
        6: 'e/bohr^3',  # Spin density
        7: 'e/bohr^3',  # WFN contribution to charge density, assuming collinear spins
        8: '1',  # Electron localization function, dimensionless
        9: 'e/bohr^3',  # Charge density minus superposition of atomic densities
        10: 'states/bohr^3',  # Integrated local density of states (ILDOS)
        11: 'Ry',  # Bare + Hartree potential
        12: 'Ry',  # the sawtooth electric field potential
        13: 'mu_B',  # Noncollinear magnetisation, Bohr magnetons
        17: 'e/bohr^3',  # All electron charge density
        18: 'T',  # The exchange and correlation magnetic field in the noncollinear case
        19: '1',  # Reduced density gradient - see dx.doi.org/10.1021/ct100641a, Eq.1 - dimensionless
        20: 'e/bohr^5',  # Product of the electron density and the second eigenvalue of the electron-density Hessian matrix, see: dx.doi.org/10.1021/ct100641a, with sign of second eigenvalue
        21: 'e/bohr^3',  # All electron charge density, PAW case
        22: 'Ry/bohr^3',  # Kinetic energy density
    }

    def parse(self, **kwargs):
        """
        Parse raw files retrieved from remote dir
        """

        temp_folder_path = None

        # A retrieved folded is required
        try:
            self.retrieved
        except exceptions.NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # If temporary files were specified, check that we have them
        if self.node.get_attribute('retrieve_temporary_list', None):
            try:
                temp_folder_path = kwargs['retrieved_temporary_folder']
            except KeyError:
                return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_TEMPORARY_FOLDER)

        # The stdout is required for parsing
        filename_stdout = self.node.get_attribute('output_filename')
        if filename_stdout not in self.retrieved.list_object_names():
            return self.exit_codes.ERROR_OUTPUT_STDOUT_MISSING

        try:
            stdout_raw = self.retrieved.get_object_content(filename_stdout)
        except (IOError, OSError):
            return self.exit_codes.ERROR_OUTPUT_STDOUT_READ

        # The post-processed data should have been written to file, either in the retrieved or temp list
        filename_data = PpCalculation._FILEOUT
        if filename_data in self.retrieved.list_object_names():   # Retrieved list case
            try:
                data_raw = self.retrieved.get_object_content(filename_data)
            except (IOError, OSError):
                return self.exit_codes.ERROR_DATAFILE_READ
        elif temp_folder_path is not None:  # Temp list case
            data_file_path = os.path.join(temp_folder_path, filename_data)
            if os.path.isfile(data_file_path):
                try:
                    with open(data_file_path, 'r') as fhandle:
                        data_raw = fhandle.read()
                except (IOError, OSError):
                    return self.exit_codes.ERROR_DATAFILE_READ
            else:
                return self.exit_codes.ERROR_OUTPUT_DATAFILE_MISSING
        else:
            return self.exit_codes.ERROR_OUTPUT_DATAFILE_MISSING

        # Parse stdout
        try:
            logs, self.output_parameters = self.parse_stdout(stdout_raw)
        except Exception:
            self.logger.error(traceback.format_exc())
            return self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION

        # Print the logs
        self.emit_logs(logs)

        # Scan logs for known errors
        if 'ERROR_PARENT_XML_MISSING' in logs['error']:
            return self.exit_codes.ERROR_PARENT_XML_MISSING
        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE' in logs['error']:
            return self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE

        # Parse the post-processed-data according to what kind of data file was produced
        if self.output_parameters['output_format'] == 'gnuplot':
            if self.output_parameters['plot_type'] == '2D polar on a sphere':
                parsed_data = self.parse_gnuplot_polar(data_raw)
            else:
                parsed_data = self.parse_gnuplot1D(data_raw)
        elif self.output_parameters['output_format'] == 'gnuplot x,y,f':
            parsed_data = self.parse_gnuplot2D(data_raw)
        elif self.output_parameters['output_format'] == 'Gaussian cube':
            parsed_data = self.parse_gaussian(data_raw)
        else:
            return self.exit_codes.ERROR_UNSUPPORTED_DATAFILE_FORMAT

        # Create output nodes
        self.out('output_data', parsed_data)
        self.out('output_parameters', orm.Dict(dict=self.output_parameters))

    def parse_stdout(self, stdout_str):
        """
        Parses the output written to StdOut to retrieve basic information about the post processing

        :param stdout_str: the stdout file read in as a single string
        """

        def detect_important_message(logs, line):
            """
            Detect know errors and warnings printed in the stdout

            :param logs:
            :param line: a line from the stdout as a string
            """
            message_map = {
                'error': {
                    'xml data file not found': 'ERROR_PARENT_XML_MISSING'
                },
                'warning': {
                    'Warning:': None,
                    'DEPRECATED:': None,
                }
            }

            # Match any known error and warning messages
            for marker, message in message_map['error'].items():
                if marker in line:
                    if message is None:
                        message = line
                    logs.error.append(message)

            for marker, message in message_map['warning'].items():
                if marker in line:
                    if message is None:
                        message = line
                    logs.warning.append(message)

        stdout_lines = stdout_str.splitlines()
        logs = get_logging_container()
        output_dict = {}

        # Check for job completion, indicating that pp.x exited without interruption, even if there was an error.
        for line in stdout_lines:
            if 'JOB DONE' in line:
                break
        else:
            logs.error.append('ERROR_OUTPUT_STDOUT_INCOMPLETE')

        # Detect any issues and detect job completion
        for line in stdout_lines:
            detect_important_message(logs, line)

        # Parse useful data from stdout
        for line in stdout_lines:
            if 'Check:' in line:
                split_line = line.split('=')
                if 'negative/imaginary' in line:    # QE6.1
                    output_dict['negative_core_charge'] = float(split_line[-1].split()[0])
                    output_dict['imaginary_core_charge'] = float(split_line.split()[-1])
                else:                               # QE6.4
                    output_dict['negative_core_charge'] = float(split_line[1])
            if 'Min, Max, imaginary charge:' in line:
                split_line = line.split()
                output_dict['charge_min'] = float(split_line[-3])
                output_dict['charge_max'] = float(split_line[-2])
                output_dict['charge_img'] = float(split_line[-1])
            if 'plot_num = ' in line:
                output_dict['plot_num'] = int(line.split('=')[1])
            if 'Plot Type:' in line:
                output_dict['plot_type'] = line.split('Output format')[0].split(':')[-1].strip()
                output_dict['output_format'] = line.split(':')[-1].strip()

        return logs, output_dict

    def parse_gnuplot1D(self, data_file_str):
        """
        Parse 1D GNUPlot formatted output

        :param data_file_str: the data file read in as a single string
        """
        data_lines = data_file_str.splitlines()

        n_col = len(data_lines[0].split())

        # 1D case
        if n_col == 2:
            coords = []
            data = []
            data_integral = []
            for line in data_lines:
                split_line = line.split()
                coords.append(float(split_line[0]))
                data.append(float(split_line[1]))
            y_data = [data]
            y_names = ['data']
            y_units = [self.units_dict[self.output_parameters['plot_num']]]

        # 1D case with spherical averaging
        if n_col == 3:
            coords = []
            data = []
            data_integral = []
            for line in data_lines:
                split_line = line.split()
                coords.append(float(split_line[0]))
                data.append(float(split_line[1]))
                data_integral.append(float(split_line[2]))
            y_data = [data, data_integral]
            y_names = ['data', 'integrated_data']
            unit = self.units_dict[self.output_parameters['plot_num']]
            y_units = [unit, unit.replace('bohr^3', 'bohr')]

        x_units = 'bohr'
        arraydata = orm.ArrayData()
        arraydata.set_array('x_coordinates', np.array(coords))
        arraydata.set_array('x_coordinates_units', np.array(x_units))
        for name, data, units in zip(y_names, y_data, y_units):
            arraydata.set_array(name, np.array(data))
            arraydata.set_array(name + '_units', np.array(units))

        return arraydata

    def parse_gnuplot_polar(self, data_file_str):
            """
            Parse 2D Polar GNUPlot formatted, single column output

            :param data_file_str: the data file read in as a single string
            """
            data_lines = data_file_str.splitlines()
            data_lines.pop(0)  # First line is a header

            data = []
            for line in data_lines:
                data.append(float(line))
            data_units = [self.units_dict[self.output_parameters['plot_num']]]

            arraydata = orm.ArrayData()
            arraydata.set_array('data', np.array(data))
            arraydata.set_array('data_units', np.array(data_units))

            return arraydata

    def parse_gnuplot2D(self, data_file_str):
        """
        Parse 2D GNUPlot formatted output

        :param data_file_str: the data file read in as a single string
        """
        data_lines = data_file_str.splitlines()

        coords = []
        data = []

        for line in data_lines:
            if line == '':
                continue
            else:
                split_line = line.split()
                coords.append([float(split_line[0]), float(split_line[1])])
                data.append(float(split_line[2]))

        coords_units = 'bohr'
        data_units = self.units_dict[self.output_parameters['plot_num']]
        arraydata = orm.ArrayData()
        arraydata.set_array('xy_coordinates', np.array(coords))
        arraydata.set_array('data', np.array(data))
        arraydata.set_array('xy_coordinates_units', np.array(coords_units))
        arraydata.set_array('data_units', np.array(data_units))

        return arraydata

    def parse_gaussian(self, data_file_str):
        """
        Parse Gaussian Cube formatted output

        :param data_file_str: the data file read in as a single string
        """

        lines = data_file_str.splitlines()

        atoms_line = lines[2].split()
        atoms = int(atoms_line[0])  # The number of atoms listed in the file
        header = lines[:6 + atoms]  # The header of the file: comments, the voxel, and the number of atoms and datapoints
        data_lines = lines[6 + atoms:]  # The actual data: atoms and volumetric data

        # Parse the declared dimensions of the volumetric data
        x_line = header[3].split()
        xdim = int(x_line[0])
        y_line = header[4].split()
        ydim = int(y_line[0])
        z_line = header[5].split()
        zdim = int(z_line[0])

        # Get the vectors describing the basis voxel
        voxel_array = np.array(
            [[x_line[1], x_line[2], x_line[3]],
            [y_line[1], y_line[2], y_line[3]],
            [z_line[1], z_line[2], z_line[3]]],
            dtype=np.float64
        )

        # Get the volumetric data
        data_array = np.zeros((xdim, ydim, zdim))

        # The data is organised into columns with a fixed number of values
        # Gather up all the data points into one list for easier unpacking
        datalist = []
        for line in data_lines:
            for i in range(0, len(line), 13):
                data_point = line[i:i + 13].strip()
                if data_point != '':
                    datalist.append(float(data_point))

        # Unpack the list and repack as a 3D array
        # Note the unusual indexing: cube files run over the z index first, then y and x.
        # E.g. The first volumetric data point is x,y,z = (0,0,0) and the second is (0,0,1)
        for i in range(0, xdim):
            for j in range(0, ydim):
                for k in range(0, zdim):
                    data_array[i, j, k] = (datalist[(i * ydim * zdim) + (j * zdim) + k])

        coordinates_units = 'bohr'
        data_units = self.units_dict[self.output_parameters['plot_num']]

        arraydata = orm.ArrayData()
        arraydata.set_array('voxel', voxel_array)
        arraydata.set_array('data', data_array)
        arraydata.set_array('coordinates_units', np.array(coordinates_units))
        arraydata.set_array('data_units', np.array(data_units))

        return arraydata
