# -*- coding: utf-8 -*-
from __future__ import absolute_import
import numpy as np
from aiida.parsers import Parser
from aiida.orm import Dict, XyData
from aiida.common import NotExistent
from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers import parse_raw_out_basic
from six.moves import range


class DosParser(Parser):
    """
    This class is the implementation of the Parser class for Dos.
    """
    # deprecated:
    #_dos_name = 'output_dos'
    #_units_name = 'output_units'


    def parse(self, **kwargs):
        """
        Parses the datafolder, stores results.
        Retrieves dos output, and some basic information from the
        out_file, such as warnings and wall_time
        """
        try:
            out_folder = self.retrieved
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # Read standard out
        try:
            filename_stdout = self.node.get_option('output_filename')  # or get_attribute(), but this is clearer
            with out_folder.open(filename_stdout, 'r') as fil:
                    out_file = fil.readlines()
        except OSError:
            return self.exit_codes.ERROR_READING_OUTPUT_FILE

        job_done = False
        for i in range(len(out_file)):
            line = out_file[-i]
            if "JOB DONE" in line:
                job_done = True
                break
        if not job_done:
            return self.exit_codes.ERROR_JOB_NOT_DONE

        # check that the dos file is present, if it is, read it
        try:
            with out_folder.open(self.node.process_class._DOS_FILENAME, 'r') as fil:
                    dos_file = fil.readlines()
        except OSError:
            return self.exit_codes.ERROR_READING_DOS_FILE

        # end of initial checks

        array_names = [[], []]
        array_units = [[], []]
        array_names[0] = ['dos_energy', 'dos',
                          'integrated_dos']  # When spin is not displayed
        array_names[1] = ['dos_energy', 'dos_spin_up', 'dos_spin_down',
                          'integrated_dos']  # When spin is displayed
        array_units[0] = ['eV', 'states/eV',
                          'states']  # When spin is not displayed
        array_units[1] = ['eV', 'states/eV', 'states/eV',
                          'states']  # When spin is displayed

        # grabs parsed data from aiida.dos
        # TODO: should I catch any QEOutputParsingError from parse_raw_dos,
        #       log an error and return an exit code?
        array_data, spin = parse_raw_dos(dos_file, array_names, array_units)
        
        energy_units = 'eV'
        dos_units = 'states/eV'
        int_dos_units = 'states'        
        xy_data = XyData()
        xy_data.set_x(array_data["dos_energy"],"dos_energy", energy_units)
        y_arrays = []
        y_names = []
        y_units = []
        y_arrays  += [array_data["integrated_dos"]]
        y_names += ["integrated_dos"]
        y_units += ["states"]
        if spin:
            y_arrays  += [array_data["dos_spin_up"]]
            y_arrays  += [array_data["dos_spin_down"]]
            y_names += ["dos_spin_up"]
            y_names += ["dos_spin_down"]
            y_units += ["states/eV"]*2
        else:
            y_arrays  += [array_data["dos"]]
            y_names += ["dos"]
            y_units += ["states/eV"]
        xy_data.set_y(y_arrays,y_names,y_units)

        # grabs the parsed data from aiida.out
        parsed_data = parse_raw_out_basic(out_file, "DOS")
        output_params = Dict(dict=parsed_data)
        # Adds warnings
        for message in parsed_data['warnings']:
            self.logger.error(message)
        # Create New Nodes
        self.out('output_dos', xy_data)
        self.out('output_parameters', output_params)


def parse_raw_dos(dos_file, array_names, array_units):
    """
    This function takes as input the dos_file as a list of filelines along
    with information on how to give labels and units to the parsed data
    
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
        raise QEOutputParsingError('dosfile could not be loaded '
        ' using genfromtxt')
    if len(dos_data) == 0:
        raise QEOutputParsingError("Dos file is empty.")
    if np.isnan(dos_data).any():
        raise QEOutputParsingError("Dos file contains non-numeric elements.")

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
        raise QEOutputParsingError("Dos file in format that the parser is not "
                                   "designed to handle.")

    i = 0
    array_data = {}
    array_data['header'] = np.array(dos_header)
    while i < len(array_names):
        array_data[array_names[i]] = dos_data[:, i]
        array_data[array_names[i]+'_units'] = np.array(array_units[i])
        i += 1
    return array_data, spin
