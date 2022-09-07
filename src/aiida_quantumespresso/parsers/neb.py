# -*- coding: utf-8 -*-
from aiida.common import NotExistent
from aiida.orm import Dict

from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.parse_raw import convert_qe_to_aiida_structure
from aiida_quantumespresso.parsers.parse_raw.neb import parse_raw_output_neb
from aiida_quantumespresso.parsers.parse_raw.pw import parse_stdout as parse_pw_stdout
from aiida_quantumespresso.parsers.parse_raw.pw import reduce_symmetries
from aiida_quantumespresso.parsers.parse_xml.exceptions import XMLParseError, XMLUnsupportedFormatError
from aiida_quantumespresso.parsers.parse_xml.pw.parse import parse_xml as parse_pw_xml
from aiida_quantumespresso.parsers.pw import PwParser

from .base import Parser


class NebParser(Parser):
    """`Parser` implementation for the `NebCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed `NebCalculation` into output nodes.

        Two nodes that are expected are the default 'retrieved' `FolderData` node which will store the retrieved files
        permanently in the repository. The second required node is a filepath under the key `retrieved_temporary_files`
        which should contain the temporary retrieved files.
        """
        import os

        from aiida.orm import ArrayData, TrajectoryData
        import numpy

        PREFIX = self.node.process_class._PREFIX

        retrieved = self.retrieved
        list_of_files = retrieved.base.repository.list_object_names()  # Note: this includes folders, but not the files they contain.

        # The stdout is required for parsing
        filename_stdout = self.node.base.attributes.get('output_filename')

        if filename_stdout not in list_of_files:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_READ)

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
        pw_input_dict = self.node.inputs.pw.parameters.get_dict()

        # load the neb input parameters dictionary
        neb_input_dict = self.node.inputs.parameters.get_dict()

        # First parse the Neb output
        try:
            stdout = retrieved.base.repository.get_object_content(filename_stdout)
            neb_out_dict, iteration_data, raw_successful = parse_raw_output_neb(stdout, neb_input_dict)
            # TODO: why do we ignore raw_successful ?
        except (OSError, QEOutputParsingError):
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_READ)

        for warn_type in ['warnings', 'parser_warnings']:
            for message in neb_out_dict[warn_type]:
                self.logger.warning(f'parsing NEB output: {message}')

        if 'QE neb run did not reach the end of the execution.' in neb_out_dict['parser_warnings']:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE)

        # Retrieve the number of images
        try:
            num_images = neb_input_dict['num_of_images']
        except KeyError:
            try:
                num_images = neb_out_dict['num_of_images']
            except KeyError:
                return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_PARSE)
        if num_images < 2:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_PARSE)

        # Now parse the information from the individual pw calculations for the different images
        image_data = {}
        positions = []
        cells = []
        # for each image...
        for i in range(num_images):
            # check if any of the known XML output file names are present, and parse the first that we find
            relative_output_folder = os.path.join(f'{PREFIX}_{i + 1}', f'{PREFIX}.save')
            retrieved_files = self.retrieved.base.repository.list_object_names(relative_output_folder)
            for xml_filename in PwCalculation.xml_filenames:
                if xml_filename in retrieved_files:
                    xml_file_path = os.path.join(relative_output_folder, xml_filename)
                    try:
                        with retrieved.base.repository.open(xml_file_path) as xml_file:
                            parsed_data_xml, logs_xml = parse_pw_xml(xml_file, None)
                    except IOError:
                        return self.exit(self.exit_codes.ERROR_OUTPUT_XML_READ)
                    except XMLParseError:
                        return self.exit(self.exit_codes.ERROR_OUTPUT_XML_PARSE)
                    except XMLUnsupportedFormatError:
                        return self.exit(self.exit_codes.ERROR_OUTPUT_XML_FORMAT)
                    except Exception as exc:
                        import traceback
                        traceback.print_exc()
                        return self.exit(self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION.format(exception=exc))
                    # this image is dealt with, so break the inner loop and go to the next image
                    break
            # otherwise, if none of the filenames we tried exists, exit with an error
            else:
                return self.exit(self.exit_codes.ERROR_MISSING_XML_FILE)

            # look for pw output and parse it
            pw_out_file = os.path.join(f'{PREFIX}_{i + 1}', 'PW.out')
            try:
                with retrieved.base.repository.open(pw_out_file, 'r') as f:
                    pw_out_text = f.read()  # Note: read() and not readlines()
            except IOError:
                return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_READ)

            try:
                parsed_data_stdout, logs_stdout = parse_pw_stdout(
                    pw_out_text, pw_input_dict, parser_options, parsed_data_xml
                )
            except Exception as exc:
                return self.exit(self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION.format(exception=exc))

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

            structure_data = convert_qe_to_aiida_structure(parsed_structure)

            key = f'pw_output_image_{i + 1}'
            image_data[key] = parsed_parameters

            positions.append([site.position for site in structure_data.sites])
            cells.append(structure_data.cell)

            # Add also PW warnings and errors to the neb output data, avoiding repetitions.
            for log_type in ['warning', 'error']:
                for message in logs_stdout[log_type]:
                    formatted_message = f'{log_type}: {message}'
                    if formatted_message not in neb_out_dict['warnings']:
                        neb_out_dict['warnings'].append(formatted_message)

        # Symbols can be obtained simply from the last image
        symbols = [str(site.kind_name) for site in structure_data.sites]

        output_params = Dict(dict(list(neb_out_dict.items()) + list(image_data.items())))
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
                for k, v in iteration_data.items():
                    arraydata.set_array(k, numpy.array(v))
                self.out('iteration_array', arraydata)

        # Load the original and interpolated energy profile along the minimum-energy path (mep)
        try:
            filename = PREFIX + '.dat'
            with retrieved.base.repository.open(filename, 'r') as handle:
                mep = numpy.loadtxt(handle)
        except Exception:
            self.logger.warning(f'could not open expected output file `{filename}`.')
            mep = numpy.array([[]])

        try:
            filename = PREFIX + '.int'
            with retrieved.base.repository.open(filename, 'r') as handle:
                interp_mep = numpy.loadtxt(handle)
        except Exception:
            self.logger.warning(f'could not open expected output file `{filename}`.')
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
