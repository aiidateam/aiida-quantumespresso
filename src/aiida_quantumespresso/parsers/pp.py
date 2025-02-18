# -*- coding: utf-8 -*-
"""`Parser` implementation for the `PpCalculation` calculation job class."""
import os
import re
from typing import Tuple

from aiida import orm
from aiida.common import AttributeDict
import numpy as np

from aiida_quantumespresso.calculations.pp import PpCalculation
from aiida_quantumespresso.utils.mapping import get_logging_container

from .base import BaseParser


class PpParser(BaseParser):
    """``Parser`` implementation for the ``PpCalculation`` calculation job class."""

    class_error_map = {
        'xml data file not found': 'ERROR_PARENT_XML_MISSING',
    }

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
        20:
        'e/bohr^5',  # Product of the electron density and the second eigenvalue of the electron-density Hessian matrix, see: dx.doi.org/10.1021/ct100641a, with sign of second eigenvalue
        21: 'e/bohr^3',  # All electron charge density, PAW case
        22: 'Ry/bohr^3',  # Kinetic energy density
    }

    def parse(self, **kwargs):
        """Parse the retrieved files of a ``PpCalculation`` into output nodes."""
        logs = get_logging_container()

        stdout, parsed_data, logs = self.parse_stdout_from_retrieved(logs)

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        parsed_pp, logs = self.parse_stdout(stdout, logs)
        parsed_data.update(parsed_pp)

        self.out('output_parameters', orm.Dict(parsed_data))

        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE'in logs.error:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE, logs)

        retrieve_temporary_list = self.node.base.attributes.get('retrieve_temporary_list', None)

        # If temporary files were specified, check that we have them
        if retrieve_temporary_list:
            try:
                retrieved_temporary_folder = kwargs['retrieved_temporary_folder']
            except KeyError:
                return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_TEMPORARY_FOLDER)

        # Currently all plot output files should start with the `filplot` as prefix. If only one file was produced the
        # prefix is the entire filename, but in the case of multiple files, there will be pairs of two files where the
        # first has the format '{filename_prefix}.{some_random_suffix' and the second has the same name but with the
        # `filename_suffix` appended.
        filename_prefix = PpCalculation._FILPLOT
        filename_suffix = PpCalculation._FILEOUT

        # How to get the output filenames and how to open them, depends on whether they will have been retrieved in the
        # `retrieved` output node, or in the `retrieved_temporary_folder`. Instead of having a conditional with almost
        # the same loop logic in each branch, we apply a somewhat dirty trick to define an `opener` which is a callable
        # that will open a handle to the output file given a certain filename. This works since it is guaranteed that
        # these output files (excluding the standard output) will all either be in the retrieved, or in the retrieved
        # temporary folder.
        if retrieve_temporary_list:
            filenames = os.listdir(retrieved_temporary_folder)
            file_opener = lambda filename: open(os.path.join(retrieved_temporary_folder, filename))
        else:
            filenames = self.retrieved.base.repository.list_object_names()
            file_opener = self.retrieved.base.repository.open

        # The following check should in principle always succeed since the iflag should in principle be set by the
        # `PpCalculation` plugin which only ever sets 0 - 4, but we check in order for the code not to except.
        iflag = self.node.inputs.parameters.base.attributes.get('PLOT')['iflag']
        if iflag not in range(5):
            return self.exit_codes.ERROR_UNSUPPORTED_DATAFILE_FORMAT

        data_parsed = []
        parsers = {
            0: self.parse_gnuplot1D,
            1: self.parse_gnuplot1D,
            2: self.parse_gnuplot2D,
            3: self.parse_gaussian,
            4: self.parse_gnuplot_polar,
        }

        def get_key_from_filename(filename):
            """Determine the output link label for the output file with the given filename."""
            if filename == filename_suffix:
                return filename

            # Pattern to match both "aiida.filplot_XXXaiida.fileout" and "aiida.filplotXXXaiida.fileout"
            pattern = rf'{re.escape(filename_prefix)}_?(.*?){re.escape(filename_suffix)}'
            matches = re.search(pattern, filename)

            if matches:
                return matches.group(1).rstrip('_')


        if self.node.base.attributes.get('parse_data_files'):
            for filename in filenames:
                # Directly parse the retrieved files after reading them to memory (`data_raw`). The raw data
                # of each file is released from memory after parsing, to improve memory usage.
                if filename.endswith(filename_suffix):
                    # Read the file to memory
                    try:
                        with file_opener(filename) as handle:
                            data_raw = handle.read()
                    except OSError:
                        return self.exit_codes.ERROR_OUTPUT_DATAFILE_READ.format(filename=filename)
                    # Parse the file
                    try:
                        key = get_key_from_filename(filename)
                        data_parsed.append((key, parsers[iflag](data_raw, self.units_dict[parsed_data['plot_num']])))
                        del data_raw
                    except Exception as exception:  # pylint: disable=broad-except
                        return self.exit_codes.ERROR_OUTPUT_DATAFILE_PARSE.format(filename=filename, exception=exception)

            # If we don't have any parsed files, we exit. Note that this will not catch the case where there should be more
            # than one file, but the engine did not retrieve all of them. Since often we anyway don't know how many files
            # should be retrieved there really is no way to check this explicitly.
            if not data_parsed:
                return self.exit_codes.ERROR_OUTPUT_DATAFILE_MISSING.format(filename=filename_prefix)

            if len(data_parsed) == 1:
                self.out('output_data', data_parsed[0][1])
            else:
                self.out('output_data_multiple', dict(data_parsed))

        return self.exit(logs=logs)

    def parse_stdout(self, stdout: str, logs: AttributeDict) -> Tuple[dict, AttributeDict]:
        """Parse the ``stdout`` content of a Quantum ESPRESSO ``pp.x`` calculation."""
        parsed_data = {}

        # Parse useful data from stdout
        for line in stdout.splitlines():
            if 'Check:' in line:  # QE < 6.5
                split_line = line.split('=')
                if 'negative/imaginary' in line:  # QE6.1-6.3
                    parsed_data['negative_core_charge'] = float(split_line[-1].split()[0])
                    parsed_data['imaginary_core_charge'] = float(split_line[-1].split()[-1])
                else:  # QE6.4
                    parsed_data['negative_core_charge'] = float(split_line[1])
            if 'Min, Max, imaginary charge:' in line:
                split_line = line.split()
                parsed_data['charge_min'] = float(split_line[-3])
                parsed_data['charge_max'] = float(split_line[-2])
                parsed_data['charge_img'] = float(split_line[-1])
            if 'plot_num = ' in line:
                parsed_data['plot_num'] = int(line.split('=')[1])
            if 'Plot Type:' in line:
                parsed_data['plot_type'] = line.split('Output format')[0].split(':')[-1].strip()
                parsed_data['output_format'] = line.split(':')[-1].strip()

        return parsed_data, logs

    @staticmethod
    def parse_gnuplot1D(data_file_str, data_units):
        """Parse 1D GNUPlot formatted output.

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
            y_units = [data_units]

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
            y_units = [data_units, data_units.replace('bohr^3', 'bohr')]

        x_units = 'bohr'
        arraydata = orm.ArrayData()
        arraydata.set_array('x_coordinates', np.array(coords))
        arraydata.set_array('x_coordinates_units', np.array(x_units))
        for name, data, units in zip(y_names, y_data, y_units):
            arraydata.set_array(name, np.array(data))
            arraydata.set_array(name + '_units', np.array(units))

        return arraydata

    @staticmethod
    def parse_gnuplot_polar(data_file_str, data_units):
        """Parse 2D Polar GNUPlot formatted, single column output.

        :param data_file_str: the data file read in as a single string
        """
        data_lines = data_file_str.splitlines()
        data_lines.pop(0)  # First line is a header

        data = []
        for line in data_lines:
            data.append(float(line))

        arraydata = orm.ArrayData()
        arraydata.set_array('data', np.array(data))
        arraydata.set_array('data_units', np.array([data_units]))

        return arraydata

    @staticmethod
    def parse_gnuplot2D(data_file_str, data_units):
        """Parse 2D GNUPlot formatted output.

        :param data_file_str: the data file read in as a single string
        """
        data_lines = data_file_str.splitlines()

        coords = []
        data = []

        for line in data_lines:
            stripped = line.strip()
            if stripped == '':
                continue
            else:
                split_line = stripped.split()
                coords.append([float(split_line[0]), float(split_line[1])])
                data.append(float(split_line[2]))

        coords_units = 'bohr'
        arraydata = orm.ArrayData()
        arraydata.set_array('xy_coordinates', np.array(coords))
        arraydata.set_array('data', np.array(data))
        arraydata.set_array('xy_coordinates_units', np.array(coords_units))
        arraydata.set_array('data_units', np.array(data_units))

        return arraydata

    @staticmethod
    def parse_gaussian(data_file_str, data_units):
        """Parse Gaussian Cube formatted output.

        :param data_file_str: the data file read in as a single string
        """
        lines = data_file_str.splitlines()

        atoms_line = lines[2].split()
        natoms = int(atoms_line[0])  # The number of atoms listed in the file
        origin = np.array(atoms_line[1:], dtype=float)

        header = lines[:6 + natoms]  # Header of the file: comments, the voxel, and the number of atoms and datapoints
        data_lines = lines[6 + natoms:]  # The actual data: atoms and volumetric data

        # Parse the declared dimensions of the volumetric data
        x_line = header[3].split()
        xdim = int(x_line[0])
        y_line = header[4].split()
        ydim = int(y_line[0])
        z_line = header[5].split()
        zdim = int(z_line[0])

        # Get the vectors describing the basis voxel
        voxel_array = np.array([[x_line[1], x_line[2], x_line[3]], [y_line[1], y_line[2], y_line[3]],
                                [z_line[1], z_line[2], z_line[3]]],
                               dtype=np.float64)

        # Get the volumetric data
        data_array = np.empty(xdim * ydim * zdim, dtype=float)
        cursor = 0
        for line in data_lines:
            ls = line.split()
            data_array[cursor:cursor + len(ls)] = ls
            cursor += len(ls)
        data_array = data_array.reshape((xdim, ydim, zdim))

        coordinates_units = 'bohr'

        arraydata = orm.ArrayData()
        arraydata.set_array('voxel', voxel_array)
        arraydata.set_array('data', data_array)
        arraydata.set_array('data_units', np.array(data_units))
        arraydata.set_array('coordinates_units', np.array(coordinates_units))

        return arraydata
