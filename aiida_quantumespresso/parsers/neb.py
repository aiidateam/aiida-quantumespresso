# -*- coding: utf-8 -*-
from __future__ import absolute_import
from aiida.parsers.parser import Parser
from aiida.orm.nodes.data.array.bands import KpointsData
from aiida.orm.nodes.data.dict import Dict
from aiida_quantumespresso.parsers import convert_qe2aiida_structure
from aiida_quantumespresso.parsers.parse_xml.pw.parse import parse_pw_xml_post_6_2
from aiida_quantumespresso.parsers.parse_xml.pw.versions import get_xml_file_version, QeXmlVersion
from aiida_quantumespresso.parsers.parse_raw.pw import parse_pw_xml_pre_6_2, parse_pw_text_output
from aiida_quantumespresso.parsers.parse_raw.pw import QEOutputParsingError
from aiida_quantumespresso.parsers.parse_raw.neb import parse_raw_output_neb
from aiida_quantumespresso.calculations.neb import NebCalculation
import six
from six.moves import range


class NebParser(Parser):
    """
    This class is the implementation of the Parser class for Neb.
    """

    def parse(self, **kwargs):
        """
        Parses the calculation-output datafolder, and stores
        results.
        """
        from aiida.common import InvalidOperation
        from aiida.orm import TrajectoryData, ArrayData
        import os
        import numpy
        import copy

        PREFIX = self.node.process_class._PREFIX

        # Check that the retrieved folder is there
        try:
            out_folder = self.retrieved
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        list_of_files = out_folder.list_object_names()  # Note: this includes folders, but not the files they contain.

        # The stdout is required for parsing
        filename_stdout = self.node.get_attribute('output_filename')

        if filename_stdout not in list_of_files:
            self.logger.error("The standard output file '{}' was not found but is required".format(filename_stdout))
            return self.exit_codes.ERROR_READING_OUTPUT_FILE

        # Look for optional settings input node and potential 'parser_options' dictionary within it
        # Note that we look for both NEB and PW parser options under "inputs.settings.parser_options";
        # we don't even have a namespace "inputs.pw.settings".
        try:
            settings = self.node.inputs.settings.get_dict()
            parser_options = settings[self.get_parser_settings_key()]
        except (AttributeError, KeyError, NotExistent):
            settings = {}
            parser_options = {}

        # load the pw input parameters dictionary
        pw_input_dict = self.node.inputs.pw__parameters.get_dict()

        # load the neb input parameters dictionary
        neb_input_dict = self.node.inputs.parameters.get_dict()

        filepath_stdout = os.path.join(out_folder._repository._get_base_folder().abspath, filename_stdout)

        # First parse the Neb output
        try:
            neb_out_dict, iteration_data, raw_successful = parse_raw_output_neb(filepath_stdout, neb_input_dict)
        except QEOutputParsingError as exc:
            self.logger.error("QEOutputParsingError in parse_raw_output_neb: {}".format(exc))
            return self.exit_codes.ERROR_READING_OUTPUT_FILE

        # Retrieve the number of images
        try:
            num_images = neb_input_dict['num_of_images']
        except KeyError:
            try:
                num_images = neb_out_dict['num_of_images']
            except KeyError:
                self.logger.error("Impossible to understand the number of images")
                return self.exit_codes.ERROR_INVALID_OUTPUT
        if num_images < 2:
            self.logger.error("Too few images: {}".format(num_images))
            return self.exit_codes.ERROR_INVALID_OUTPUT
        
        # Now parse the information from the single pw calculations for the different images
        image_data = {}
        positions = []
        cells = []
        # for each image...
        for i in range(num_images):
            # check if any of the known XML output file names are present, and parse the first that we find
            relative_output_folder = os.path.join('{}_{}'.format(PREFIX, i+1), '{}.save'.format(PREFIX))
            retrieved_files = self.retrieved.list_object_names(relative_output_folder)
            for xml_filename in PwCalculation.xml_filenames:
                if xml_filename in retrieved_files:
                    xml_file_abspath = os.path.join(out_folder._repository._get_base_folder().abspath,
                                                    relative_output_folder,
                                                    xml_filename)
                    try:
                        xml_file_version = get_xml_file_version(xml_file_abspath)
                    except ValueError as exception:
                        raise QEOutputParsingError('failed to determine XML output file version: {}'.format(exception))

                    if xml_file_version == QeXmlVersion.POST_6_2:
                        xml_data, structure_dict, bands_data = parse_pw_xml_post_6_2(xml_file_abspath, parser_options, self.logger)
                    elif xml_file_version == QeXmlVersion.PRE_6_2:
                        xml_data, structure_dict, bands_data = parse_pw_xml_pre_6_2(xml_file_abspath, None, parser_options, self.logger)
                    else:
                        raise ValueError('unrecognized XML file version')
                    
                    structure_data = convert_qe2aiida_structure(structure_dict)
                    break
            # otherwise, if none of the filenames we tried exists, exit with an error
            else:
                self.logger.error("No xml output file found for image {}".format(i+1))
                return self.exit_codes.ERROR_MISSING_XML_FILE
            
            # look for pw output and parse it
            pw_out_file = os.path.join(out_folder._repository._get_base_folder().abspath,
                                       '{}_{}'.format(PREFIX, i+1),
                                       'PW.out')
            try:
                with open(pw_out_file,'r') as f:
                    pw_out_lines = f.read()  # Note: read() and not readlines()
            except IOError:
                self.logger.error("No pw output file found for image {}".format(i+1))
                return self.exit_codes.ERROR_READING_OUTPUT_FILE
            
            pw_out_data, trajectory_data, critical_messages = parse_pw_text_output(pw_out_lines, xml_data,
                                                                                   structure_dict, pw_input_dict)
            
            # I add in the out_data all the last elements of trajectory_data values.
            # Safe for some large arrays, that I will likely never query.
            skip_keys = ['forces','atomic_magnetic_moments','atomic_charges',
                         'lattice_vectors_relax','atomic_positions_relax',
                         'atomic_species_name']
            tmp_trajectory_data = copy.copy(trajectory_data)
            for x in six.iteritems(tmp_trajectory_data):
                if x[0] in skip_keys:
                    continue
                pw_out_data[x[0]] = x[1][-1]
                if len(x[1]) == 1:  # delete any keys that are not arrays
                    trajectory_data.pop(x[0])
            # As the k points are an array that is rather large, and again it's not something I'm going to parse likely
            # since it's an info mainly contained in the input file, I move it to the trajectory data
            for key in ['k_points','k_points_weights']:
                try:
                    trajectory_data[key] = xml_data.pop(key)
                except KeyError:
                    pass

            key = 'pw_output_image_{}'.format(i+1)
            image_data[key] = dict(list(pw_out_data.items()) + list(xml_data.items()))

            positions.append([site.position for site in structure_data.sites])
            cells.append(structure_data.cell)

            # If a warning was already present in the NEB, add also PW warnings to the neb output data,
            # avoiding repetitions.
            if neb_out_dict['warnings']:
                for warning in pw_out_data['warnings']:
                    if warning not in neb_out_dict['warnings']:
                        neb_out_dict['warnings'].append(warning) 
        
        # Symbols can be obtained simply from the last image
        symbols = [str(site.kind_name) for site in structure_data.sites]
        
        # convert output parameters to AiiDA Dict, and add as output node
        output_params = Dict(dict=dict(list(neb_out_dict.items())+list(image_data.items())))
        self.out('output_parameters', output_params)
        
        # convert data on structure of images into a TrajectoryData, and add as output node
        traj = TrajectoryData()
        traj.set_trajectory(stepids = numpy.arange(1,num_images+1),
                            cells = numpy.array(cells),
                            symbols = numpy.array(symbols),
                            positions = numpy.array(positions),
                            )
        self.out('output_trajectory', traj)
        
        if parser_options.get('all_iterations', False):
            if iteration_data:
                arraydata = ArrayData()
                for k, v in six.iteritems(iteration_data):
                    arraydata.set_array(k, numpy.array(v))
                self.out('iteration_array', arraydata)
        
        # Load the original and interpolated energy profile along the minimum-energy path (mep)
        try:
            mep_file = os.path.join( out_folder._repository._get_base_folder().abspath,
                                     PREFIX + '.dat' )
            mep = numpy.loadtxt(mep_file)
        except Exception:
            self.logger.warning("Impossible to find the file with image energies "
                                "versus reaction coordinate.")
            mep = numpy.array([[]])

        try:
            interp_mep_file = os.path.join( out_folder._repository._get_base_folder().abspath,
                                            PREFIX + '.int' )
            interp_mep = numpy.loadtxt(interp_mep_file)
        except Exception:
            self.logger.warning("Impossible to find the file with the interpolation "
                                "of image energies versus reaction coordinate.")
            interp_mep = numpy.array([[]])
        
        # Create an ArrayData with the energy profiles
        mep_arraydata = ArrayData()
        mep_arraydata.set_array('mep', mep)
        mep_arraydata.set_array('interpolated_mep', interp_mep)
        self.out('output_mep', mep_arraydata)
        
        return

    def get_parser_settings_key(self):
        """
        Return the name of the key to be used in the calculation settings, that
        contains the dictionary with the parser_options 
        """
        return 'parser_options'
