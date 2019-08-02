# -*- coding: utf-8 -*-
from __future__ import absolute_import

import six
from six.moves import range

from aiida.parsers.parser import Parser
from aiida.orm import Dict, KpointsData
from aiida.common import NotExistent
from aiida_quantumespresso.parsers import convert_qe2aiida_structure, QEOutputParsingError
from aiida_quantumespresso.parsers.parse_raw.pw import parse_stdout as parse_pw_stdout
from aiida_quantumespresso.parsers.parse_xml.pw.parse import parse_xml as parse_pw_xml
from aiida_quantumespresso.parsers.parse_xml.pw.exceptions import XMLParseError, XMLUnsupportedFormatError
from aiida_quantumespresso.parsers.parse_raw.neb import parse_raw_output_neb
from aiida_quantumespresso.parsers.pw import PwParser
from aiida_quantumespresso.calculations.pw import PwCalculation


class NebParser(Parser):
    """
    This class is the implementation of the Parser class for Neb.
    """

    def parse(self, **kwargs):
        """
        Parses the calculation-output datafolder, and stores
        results.
        """
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

        try:
            include_deprecated_v2_keys = parser_options['include_deprecated_v2_keys']
        except (TypeError,KeyError):
            include_deprecated_v2_keys = False

        # load the pw input parameters dictionary
        pw_input_dict = self.node.inputs.pw__parameters.get_dict()

        # load the neb input parameters dictionary
        neb_input_dict = self.node.inputs.parameters.get_dict()

        filepath_stdout = os.path.join(out_folder._repository._get_base_folder().abspath, filename_stdout)

        # First parse the Neb output
        try:
            neb_out_dict, iteration_data, raw_successful = parse_raw_output_neb(filepath_stdout, neb_input_dict)
            # TODO: why do we ignore raw_successful ?
        except QEOutputParsingError as exc:
            self.logger.error('QEOutputParsingError in parse_raw_output_neb: {}'.format(exc))
            return self.exit_codes.ERROR_READING_OUTPUT_FILE

        for warn_type in ['warnings', 'parser_warnings']:
            for message in neb_out_dict[warn_type]:
                self.logger.warning('parsing NEB output: {}'.format(message))

        if 'QE neb run did not reach the end of the execution.' in neb_out_dict['parser_warnings']:
            return self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE

        # Retrieve the number of images
        try:
            num_images = neb_input_dict['num_of_images']
        except KeyError:
            try:
                num_images = neb_out_dict['num_of_images']
            except KeyError:
                self.logger.error('Impossible to understand the number of images')
                return self.exit_codes.ERROR_INVALID_OUTPUT
        if num_images < 2:
            self.logger.error('Too few images: {}'.format(num_images))
            return self.exit_codes.ERROR_INVALID_OUTPUT

        # Now parse the information from the individual pw calculations for the different images
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
                        with open(xml_file_abspath) as xml_file:
                            parsed_data_xml, logs_xml = parse_pw_xml(xml_file, None, include_deprecated_v2_keys)
                    except IOError:
                        return self.exit_codes.ERROR_OUTPUT_XML_READ
                    except XMLParseError:
                        return self.exit_codes.ERROR_OUTPUT_XML_PARSE
                    except XMLUnsupportedFormatError:
                        return self.exit_codes.ERROR_OUTPUT_XML_FORMAT
                    except Exception:
                        import traceback
                        traceback.print_exc()
                        return self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION
                    # this image is dealt with, go to next image
                    break
            # otherwise, if none of the filenames we tried exists, exit with an error
            else:
                self.logger.error('No xml output file found for image {}'.format(i+1))
                return self.exit_codes.ERROR_MISSING_XML_FILE

            # look for pw output and parse it
            pw_out_file = os.path.join(out_folder._repository._get_base_folder().abspath,
                                       '{}_{}'.format(PREFIX, i+1),
                                       'PW.out')
            try:
                with open(pw_out_file,'r') as f:
                    pw_out_text = f.read()  # Note: read() and not readlines()
            except IOError:
                self.logger.error('No pw output file found for image {}'.format(i+1))
                return self.exit_codes.ERROR_READING_OUTPUT_FILE

            try:
                parsed_data_stdout, logs_stdout = parse_pw_stdout(pw_out_text, pw_input_dict, parser_options, parsed_data_xml)
            except Exception:
                import traceback
                traceback.print_exc()
                return self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION

            parsed_structure = parsed_data_stdout.pop('structure', {})
            parsed_trajectory = parsed_data_stdout.pop('trajectory', {})
            parsed_parameters = PwParser.build_output_parameters(parsed_data_xml, parsed_data_stdout)

            # Delete bands # TODO: this is just to make pytest happy; do we want to keep them instead?
            try:
                parsed_parameters.pop('bands')
            except KeyError:
                pass

            # TODO: blacklist (old, here) or whitelist (new, from PwParser - see below)?
            # # I add in the out_data all the last elements of trajectory_data values.
            # # Safe for some large arrays, that I will likely never query.
            # skip_keys = ['forces','atomic_magnetic_moments','atomic_charges',
            #              'lattice_vectors_relax','atomic_positions_relax',
            #              'atomic_species_name']
            # tmp_trajectory_data = copy.copy(parsed_trajectory)
            # for k, v in six.iteritems(tmp_trajectory_data):
            #     if k in skip_keys:
            #         continue
            #     pw_out_data[k] = v[-1]
            #     if len(v) == 1:  # delete any keys that are not arrays
            #         trajectory_data.pop(x[0])

            # As the k points are an array that is rather large, and again it's not something I'm going to parse likely
            # since it's an info mainly contained in the input file, I move it to the trajectory data.
            # Same for bands.
            # TODO: this is actually useless, since we discard the trajectory of a single PW calculation,
            #  and these 2 keys are not preserved by PwParser.final_trajectory_frame_to_parameters
            #  Do we need them?
            for key in ['k_points','k_points_weights']:
                try:
                    parsed_trajectory[key] = parsed_parameters.pop(key)
                except KeyError:
                    pass

            # Append the last frame of some of the smaller trajectory arrays to the parameters for easy querying
            PwParser.final_trajectory_frame_to_parameters(parsed_parameters, parsed_trajectory)

            # If the parser option 'all_symmetries' is not set to True, we reduce the raw parsed symmetries to safe space
            # TODO: do we need this!?
            # self.reduce_symmetries(parsed_parameters, parsed_structure, parser_options)
            # TODO: do we need these!?
            # kpoints = PwParser.build_output_kpoints(parsed_parameters, structure)
            # trajectory = PwParser.build_output_trajectory(parsed_trajectory, structure)

            structure_data = convert_qe2aiida_structure(parsed_structure)

            key = 'pw_output_image_{}'.format(i+1)
            image_data[key] = parsed_parameters

            positions.append([site.position for site in structure_data.sites])
            cells.append(structure_data.cell)

            # Add also PW warnings and errors to the neb output data, avoiding repetitions.
            for log_type in ['warning', 'error']:
                for message in logs_stdout[log_type]:
                    formatted_message = '{}: {}'.format(log_type, message)
                    if formatted_message not in neb_out_dict['warnings']:
                        neb_out_dict['warnings'].append(formatted_message)

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
            self.logger.warning('Impossible to find the file with image energies '
                                'versus reaction coordinate.')
            mep = numpy.array([[]])

        try:
            interp_mep_file = os.path.join( out_folder._repository._get_base_folder().abspath,
                                            PREFIX + '.int' )
            interp_mep = numpy.loadtxt(interp_mep_file)
        except Exception:
            self.logger.warning('Impossible to find the file with the interpolation '
                                'of image energies versus reaction coordinate.')
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
