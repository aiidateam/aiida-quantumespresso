import os

import numpy as np
from aiida.common import NotExistent
from aiida.orm import ArrayData, Dict, TrajectoryData

from aiida_quantumespresso.calculations.neb import NebCalculation
from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.parsers.parse_raw import convert_qe_to_aiida_structure
from aiida_quantumespresso.parsers.parse_raw.neb import parse_raw_output_neb
from aiida_quantumespresso.parsers.parse_raw.pw import parse_stdout as parse_pw_stdout
from aiida_quantumespresso.parsers.parse_raw.pw import reduce_symmetries
from aiida_quantumespresso.parsers.pw import PwParser
from aiida_quantumespresso.utils.mapping import get_logging_container

from .base import BaseParser


class NebParser(BaseParser):
    """`Parser` implementation for the `NebCalculation` calculation job class."""

    # Key that contains the optional parser options in the `settings` input node.
    parser_settings_key = 'parser_options'

    class_warning_map = {
        'scf convergence NOT achieved on image': 'SCF did not converge for a given image',
        'Maximum CPU time exceeded': 'Maximum CPU time exceeded',
        # !! 'step' and not 'steps' is needed in order to be found by regex
        'reached the maximum number of step': 'Maximum number of iterations reached in the image optimization',
    }

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed ``NebCalculation`` into output nodes.

        Two nodes that are expected are the default 'retrieved' ``FolderData`` node which will store the retrieved files
        permanently in the repository. The second required node is a filepath under the key `retrieved_temporary_files`
        which should contain the temporary retrieved files.
        """
        logs = get_logging_container()

        prefix = self.node.process_class._PREFIX  # noqa: SLF001

        self.exit_code_xml = None

        # Look for optional settings input node and potential 'parser_options' dictionary within it
        # Note that we look for both NEB and PW parser options under "inputs.settings.parser_options";
        # we don't even have a namespace "inputs.pw.settings".
        try:
            settings = self.node.inputs.settings.get_dict()
            parser_options = settings[self.parser_settings_key]
        except (AttributeError, KeyError, NotExistent):
            settings = {}
            parser_options = {}

        # load the pw input parameters dictionary
        pw_input_dict = self.node.inputs.pw.parameters.get_dict()

        stdout, parsed_data, logs = self.parse_stdout_from_retrieved(logs)

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        neb_out_dict, iteration_data = parse_raw_output_neb(stdout)
        parsed_data.update(neb_out_dict)

        num_images = parsed_data['num_of_images']

        # Now parse the information from the individual pw calculations for the different images
        image_data = {}
        positions = []
        cells = []

        for i in range(num_images):
            # check if any of the known XML output file names are present, and parse it
            relative_output_folder = os.path.join(f'{prefix}_{i + 1}', f'{prefix}.save')
            parsed_data_xml, logs_xml = self.parse_xml(relative_output_folder)

            # look for pw output and parse it
            pw_out_file = os.path.join(f'{prefix}_{i + 1}', 'PW.out')
            try:
                with self.retrieved.base.repository.open(pw_out_file, 'r') as f:
                    pw_out_text = f.read()  # Note: read() and not readlines()
                # Output file can contain the output of many scf iterations, analyse only the last one
                pw_out_text = '     coordinates at iteration' + pw_out_text.split('coordinates at iteration')[-1]
            except OSError:
                logs_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_READ

            try:
                parsed_data_stdout, logs_stdout = parse_pw_stdout(
                    pw_out_text, pw_input_dict, parser_options, parsed_data_xml
                )
            except Exception as exc:
                return self.exit(self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION.format(exception=exc))

            logs_stdout['error'].remove('ERROR_OUTPUT_STDOUT_INCOMPLETE')

            # Determine issues coming from electronic structure calculations
            exit_code = self.validate_electronic(logs_stdout)
            if exit_code:
                return self.exit(exit_code)

            exit_code = self.validate_premature_exit(logs_stdout)
            if exit_code:
                return self.exit(exit_code)

            if logs_stdout.error and self.exit_code_xml:
                return self.exit(self.exit_codes.ERROR_OUTPUT_FILES)

            parsed_structure = parsed_data_stdout.pop('structure', {})
            parsed_trajectory = parsed_data_xml.pop('trajectory', {})
            parsed_parameters = parsed_data_xml
            PwParser.backwards_compatibility_parameters(parsed_parameters, parsed_data_stdout)

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
            for log_level in ['warning', 'error']:
                for message in logs_stdout[log_level]:
                    formatted_message = f'{log_level}: {message}'
                    if formatted_message not in parsed_data['warnings']:
                        parsed_data['warnings'].append(formatted_message)

        # Symbols can be obtained simply from the last image
        symbols = [str(site.kind_name) for site in structure_data.sites]

        output_params = Dict(dict(list(parsed_data.items()) + list(image_data.items())))
        self.out('output_parameters', output_params)

        trajectory = TrajectoryData()
        trajectory.set_trajectory(
            stepids=np.arange(1, num_images + 1),
            cells=np.array(cells),
            symbols=symbols,
            positions=np.array(positions),
        )
        self.out('output_trajectory', trajectory)

        if parser_options is not None and parser_options.get('all_iterations', False) and iteration_data:
            arraydata = ArrayData()
            for k, v in iteration_data.items():
                arraydata.set_array(k, np.array(v))
            self.out('iteration_array', arraydata)

        # Load the original and interpolated energy profile along the minimum-energy path (mep)
        try:
            filename = prefix + '.dat'
            with self.retrieved.base.repository.open(filename, 'r') as handle:
                mep = np.loadtxt(handle)
        except Exception:
            self.logger.warning(f'could not open expected output file `{filename}`.')
            mep = np.array([[]])

        try:
            filename = prefix + '.int'
            with self.retrieved.base.repository.open(filename, 'r') as handle:
                interp_mep = np.loadtxt(handle)
        except Exception:
            self.logger.warning(f'could not open expected output file `{filename}`.')
            interp_mep = np.array([[]])

        # Create an ArrayData with the energy profiles
        mep_arraydata = ArrayData()
        mep_arraydata.set_array('mep', mep)
        mep_arraydata.set_array('interpolated_mep', interp_mep)
        self.out('output_mep', mep_arraydata)

        if logs.error:
            # First check whether the scheduler already reported an exit code.
            if self.node.exit_status is not None:
                # The following scheduler errors should correspond to cases where we can simply restart the calculation
                # and have a chance that the calculation will succeed as the error can be transient.
                recoverable_scheduler_error = self.node.exit_status in [
                    NebCalculation.exit_codes.ERROR_SCHEDULER_OUT_OF_WALLTIME.status,
                    NebCalculation.exit_codes.ERROR_SCHEDULER_NODE_FAILURE.status,
                ]
                if recoverable_scheduler_error:
                    return NebCalculation.exit_codes.ERROR_NEB_INTERRUPTED_PARTIAL_TRAJECTORY
        elif 'Maximum number of iterations reached in the image optimization' in logs.warning:
            return NebCalculation.exit_codes.ERROR_NEB_CYCLE_EXCEEDED_NSTEP
        else:
            # Calculation completed successfully shortly after exceeding walltime but before being terminated by the
            # scheduler. In that case 'exit_status' can be reset.
            self.node.set_exit_status(None)

        return self.exit(logs=logs)

    def parse_xml(self, relative_output_folder):
        """Parse the XML output file for the specific image.

        :param relative_output_folder: relative path to the output folder of the image.
        :return: tuple of two dictionaries, first with raw parsed data and second with log messages
        """
        from aiida_quantumespresso.parsers.parse_xml.exceptions import (
            XMLParseError,
            XMLUnsupportedFormatError,
        )
        from aiida_quantumespresso.parsers.parse_xml.parse import (
            parse_xml as parse_pw_xml,
        )

        logs = get_logging_container()
        parsed_data = {}

        try:
            retrieved_files = self.retrieved.base.repository.list_object_names(relative_output_folder)
        except Exception:
            self.exit_code_xml = self.exit_codes.ERROR_OUTPUT_XML_MISSING
            return parsed_data, logs

        xml_filenames = [
            os.path.join(relative_output_folder, xml_file)
            for xml_file in PwCalculation.xml_filenames
            if xml_file in retrieved_files
        ]
        if not xml_filenames:
            if not self.node.get_option('without_xml'):
                self.exit_code_xml = self.exit_codes.ERROR_OUTPUT_XML_MISSING
            return parsed_data, logs

        if len(xml_filenames) > 1:
            self.exit_code_xml = self.exit_codes.ERROR_OUTPUT_XML_MULTIPLE
            return parsed_data, logs

        try:
            with self.retrieved.base.repository.open(xml_filenames[0]) as xml_file:
                parsed_data, logs = parse_pw_xml(xml_file)
        except OSError:
            self.exit_code_xml = self.exit_codes.ERROR_OUTPUT_XML_READ
        except XMLParseError:
            self.exit_code_xml = self.exit_codes.ERROR_OUTPUT_XML_PARSE
        except XMLUnsupportedFormatError:
            self.exit_code_xml = self.exit_codes.ERROR_OUTPUT_XML_FORMAT
        except Exception:
            import traceback

            logs.critical.append(traceback.format_exc())
            self.exit_code_xml = self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION

        return parsed_data, logs

    def validate_premature_exit(self, logs):
        """Analyze problems that will cause a pre-mature termination of the calculation, controlled or not."""

        for error_label in [
            'ERROR_OUT_OF_WALLTIME',
            'ERROR_DEXX_IS_NEGATIVE',
            'ERROR_COMPUTING_CHOLESKY',
            'ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED',
            'ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE',
            'ERROR_ZHEGVD_FAILED',
            'ERROR_QR_FAILED',
            'ERROR_EIGENVECTOR_CONVERGENCE',
            'ERROR_BROYDEN_FACTORIZATION',
        ]:
            if error_label in logs['error']:
                return self.exit_codes.get(error_label)

    def validate_electronic(self, logs):
        """Analyze problems that are specific to `electronic` scf calculations."""

        if 'ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED' in logs['error']:
            scf_must_converge = self.node.inputs.pw.parameters.base.attributes.get('ELECTRONS', {}).get(
                'scf_must_converge', True
            )
            electron_maxstep = self.node.inputs.pw.parameters.base.attributes.get('ELECTRONS', {}).get(
                'electron_maxstep', 1
            )

            if electron_maxstep == 0 or not scf_must_converge:
                return self.exit_codes.WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED

            return self.exit_codes.ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED
