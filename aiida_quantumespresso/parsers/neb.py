# -*- coding: utf-8 -*-
from __future__ import absolute_import

import six
from six.moves import range

from aiida.common import NotExistent
from aiida.orm import Dict
from aiida.parsers.parser import Parser
from aiida_quantumespresso.parsers import convert_qe2aiida_structure, QEOutputParsingError
from aiida_quantumespresso.parsers.parse_raw.pw import reduce_symmetries
from aiida_quantumespresso.parsers.parse_raw.pw import parse_stdout as parse_pw_stdout
from aiida_quantumespresso.parsers.parse_xml.pw.parse import parse_xml as parse_pw_xml
from aiida_quantumespresso.parsers.parse_xml.pw.exceptions import XMLParseError, XMLUnsupportedFormatError
from aiida_quantumespresso.parsers.parse_raw.neb import parse_raw_output_neb
from aiida_quantumespresso.parsers.pw import PwParser
from aiida_quantumespresso.calculations.pw import PwCalculation


class NebParser(Parser):
    """`Parser` implementation for the `NebCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed `NebCalculation` into output nodes.

        Two nodes that are expected are the default 'retrieved' `FolderData` node which will store the retrieved files
        permanently in the repository. The second required node is a filepath under the key `retrieved_temporary_files`
        which should contain the temporary retrieved files.
        """
        from aiida.orm import TrajectoryData, ArrayData
        import os
        import numpy

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
        except (TypeError, KeyError):
            include_deprecated_v2_keys = False

        # load the pw input parameters dictionary
        pw_input_dict = self.node.inputs.pw__parameters.get_dict()

        # load the neb input parameters dictionary
        neb_input_dict = self.node.inputs.parameters.get_dict()

        stdout_abspath = os.path.join(out_folder._repository._get_base_folder().abspath, filename_stdout)

        # First parse the Neb output
        try:
            neb_out_dict, iteration_data, raw_successful = parse_raw_output_neb(stdout_abspath, neb_input_dict)
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
            relative_output_folder = os.path.join('{}_{}'.format(PREFIX, i + 1), '{}.save'.format(PREFIX))
            retrieved_files = self.retrieved.list_object_names(relative_output_folder)
            for xml_filename in PwCalculation.xml_filenames:
                if xml_filename in retrieved_files:
                    xml_file_path = os.path.join(relative_output_folder, xml_filename)
                    try:
                        with out_folder.open(xml_file_path) as xml_file:
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
                    # this image is dealt with, so break the inner loop and go to the next image
                    break
            # otherwise, if none of the filenames we tried exists, exit with an error
            else:
                self.logger.error('No xml output file found for image {}'.format(i + 1))
                return self.exit_codes.ERROR_MISSING_XML_FILE

            # look for pw output and parse it
            pw_out_file = os.path.join('{}_{}'.format(PREFIX, i + 1), 'PW.out')
            try:
                with out_folder.open(pw_out_file, 'r') as f:
                    pw_out_text = f.read()  # Note: read() and not readlines()
            except IOError:
                self.logger.error('No pw output file found for image {}'.format(i + 1))
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

            # Explicit information about k-points does not need to be queryable so we remove it from the parameters
            parsed_parameters.pop('k_points', None)
            parsed_parameters.pop('k_points_units', None)
            parsed_parameters.pop('k_points_weights', None)

            # Delete bands # TODO: this is just to make pytest happy; do we want to keep them instead?
            parsed_parameters.pop('bands', None)

            # Append the last frame of some of the smaller trajectory arrays to the parameters for easy querying
            PwParser.final_trajectory_frame_to_parameters(parsed_parameters, parsed_trajectory)

            # If the parser option 'all_symmetries' is False, we reduce the raw parsed symmetries to save space
            all_symmetries = False if parser_options is None else parser_options.get('all_symmetries', False)
            if not all_symmetries and 'cell' in parsed_structure:
                reduce_symmetries(parsed_parameters, parsed_structure, self.logger)

            structure_data = convert_qe2aiida_structure(parsed_structure)

            key = 'pw_output_image_{}'.format(i + 1)
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

        output_params = Dict(dict=dict(list(neb_out_dict.items()) + list(image_data.items())))
        self.out('output_parameters', output_params)

        trajectory = TrajectoryData()
        trajectory.set_trajectory(
            stepids=numpy.arange(1, num_images + 1),
            cells=numpy.array(cells),
            symbols=symbols,
            positions=numpy.array(positions),
        )
        self.out('output_trajectory', trajectory)

        if parser_options is not None and parser_options.get('all_iterations', False):
            if iteration_data:
                arraydata = ArrayData()
                for k, v in six.iteritems(iteration_data):
                    arraydata.set_array(k, numpy.array(v))
                self.out('iteration_array', arraydata)

        # Load the original and interpolated energy profile along the minimum-energy path (mep)
        try:
            filename = PREFIX + '.dat'
            with out_folder.open(filename, 'r') as handle:
                mep = numpy.loadtxt(handle)
        except Exception:
            self.logger.warning('could not open expected output file `{}`.'.format(filename))
            mep = numpy.array([[]])

        try:
            filename = PREFIX + '.int'
            with out_folder.open(filename, 'r') as handle:
                interp_mep = numpy.loadtxt(handle)
        except Exception:
            self.logger.warning('could not open expected output file `{}`.'.format(filename))
            interp_mep = numpy.array([[]])

        # Create an ArrayData with the energy profiles
        mep_arraydata = ArrayData()
        mep_arraydata.set_array('mep', mep)
        mep_arraydata.set_array('interpolated_mep', interp_mep)
        self.out('output_mep', mep_arraydata)

        return

    @staticmethod
    def get_parser_settings_key():
        """Return the key that contains the optional parser options in the `settings` input node."""
        return 'parser_options'
