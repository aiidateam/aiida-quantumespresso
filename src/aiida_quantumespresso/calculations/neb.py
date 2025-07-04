# -*- coding: utf-8 -*-
"""Plugin to create a Quantum Espresso neb.x input file."""
import copy
import os
import warnings

from aiida import orm
from aiida.common import CalcInfo, CodeInfo, InputValidationError
from aiida.common.lang import classproperty
from aiida.common.warnings import AiidaDeprecationWarning

from aiida_quantumespresso.calculations import _lowercase_dict, _pop_parser_options, _uppercase_dict
from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry

from .base import CalcJob


class NebCalculation(CalcJob):
    """Nudged Elastic Band code (neb.x) of Quantum ESPRESSO distribution."""

    _PREFIX = 'aiida'

    # in restarts, will not copy but use symlinks
    _default_symlink_usage = False

    # Default input and output file names
    _DEFAULT_INPUT_FILE = 'neb.dat'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'
    _PSEUDO_SUBFOLDER = PwCalculation._PSEUDO_SUBFOLDER  # pylint: disable=protected-access
    _OUTPUT_SUBFOLDER = PwCalculation._OUTPUT_SUBFOLDER  # pylint: disable=protected-access

    # Keywords that cannot be set (for the PW input)
    _blocked_keywords = []

    _use_kpoints = True

    @classproperty
    def _internal_retrieve_list(cls):
        # pylint: disable=no-self-argument
        # I retrieve them all, even if I don't parse all of them
        _neb_ext_list = ['path', 'dat', 'int']
        return [f'{cls._PREFIX}.{ext}' for ext in _neb_ext_list]

    @classproperty
    def xml_filepaths(cls):
        """Return a list of relative filepaths of XML files."""
        # pylint: disable=no-self-argument,not-an-iterable
        filepaths = []

        for filename in PwCalculation.xml_filenames:
            filepath = os.path.join(cls._OUTPUT_SUBFOLDER, cls._PREFIX + '_*[0-9]', cls._PREFIX + '.save', filename)
            filepaths.append(filepath)

        return filepaths

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('metadata.options.input_filename', valid_type=str, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=str, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.parser_name', valid_type=str, default='quantumespresso.neb')
        spec.input('first_structure', valid_type=orm.StructureData, help='Initial structure')
        spec.input('last_structure', valid_type=orm.StructureData, help='Final structure')
        spec.input('parameters', valid_type=orm.Dict, help='NEB-specific input parameters')
        spec.input('settings', valid_type=orm.Dict, required=False,
            help='Optional parameters to affect the way the calculation job and the parsing are performed.')
        spec.input('parent_folder', valid_type=orm.RemoteData, required=False,
            help='An optional working directory of a previously completed calculation to restart from.')
        # We reuse some inputs from PwCalculation to construct the PW-specific parts of the input files
        spec.expose_inputs(PwCalculation, namespace='pw', include=('parameters', 'pseudos', 'kpoints', 'vdw_table'))
        spec.output('output_parameters', valid_type=orm.Dict,
            help='The output parameters dictionary of the NEB calculation')
        spec.output('output_trajectory', valid_type=orm.TrajectoryData)
        spec.output('iteration_array', valid_type=orm.ArrayData, required=False)
        spec.output('output_mep', valid_type=orm.ArrayData,
            help='The original and interpolated energy profiles along the minimum-energy path (mep)')
        spec.default_output_node = 'output_parameters'
        spec.exit_code(303, 'ERROR_MISSING_XML_FILE',
            message='The required XML file is not present in the retrieved folder.')
        spec.exit_code(320, 'ERROR_OUTPUT_XML_READ',
            message='The XML output file could not be read.')
        spec.exit_code(321, 'ERROR_OUTPUT_XML_PARSE',
            message='The XML output file could not be parsed.')
        spec.exit_code(322, 'ERROR_OUTPUT_XML_FORMAT',
            message='The XML output file has an unsupported format.')
        spec.exit_code(350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION',
            message='The parser raised an unexpected exception: {exception}')
        # yapf: enable

    @classmethod
    def _generate_input_files(cls, neb_parameters, settings_dict):
        """Generate the input data for the NEB part of the calculation."""
        # I put the first-level keys as uppercase (i.e., namelist and card names)
        # and the second-level keys as lowercase
        # (deeper levels are unchanged)
        input_params = _uppercase_dict(neb_parameters.get_dict(), dict_name='parameters')
        input_params = {k: _lowercase_dict(v, dict_name=k) for k, v in input_params.items()}

        # Force default values for blocked keywords. NOTE: this is different from PW/CP
        for blocked in cls._blocked_keywords:
            namelist = blocked[0].upper()
            key = blocked[1].lower()
            value = blocked[2]
            if namelist in input_params:
                if key in input_params[namelist]:
                    raise InputValidationError(
                        f"You cannot specify explicitly the '{key}' key in the '{namelist}' namelist."
                    )
            else:
                input_params[namelist] = {}
            input_params[namelist][key] = value

        # Create an empty dictionary for the compulsory namelist 'PATH' if not present
        if 'PATH' not in input_params:
            input_params['PATH'] = {}

        # In case of climbing image, we need the corresponding card
        ci_scheme = input_params['PATH'].get('ci_scheme', 'no-ci').lower()
        climbing_image_list = settings_dict.pop('CLIMBING_IMAGES', None)
        if ci_scheme == 'manual':
            manual_climbing_image = True
            if climbing_image_list is None:
                raise InputValidationError(
                    f"'ci_scheme' is {ci_scheme}, but no climbing images were specified for this calculation."
                )
            if not isinstance(climbing_image_list, list):
                raise InputValidationError('Climbing images should be provided as a list.')
            num_of_images = input_params['PATH'].get('num_of_images', 2)
            if any((i < 2 or i >= num_of_images) for i in climbing_image_list):
                raise InputValidationError(
                    'The climbing images should be in the range between the first '
                    'and the last image (excluded).'
                )
            climbing_image_card = 'CLIMBING_IMAGES\n'
            climbing_image_card += ', '.join([str(_) for _ in climbing_image_list]) + '\n'
        else:
            manual_climbing_image = False
            if climbing_image_list is not None:
                raise InputValidationError(f"Climbing images are not accepted when 'ci_scheme' is {ci_scheme}.")

        input_data = '&PATH\n'
        # namelist content; set to {} if not present, so that we leave an empty namelist
        namelist = input_params.pop('PATH', {})
        for key, value in sorted(namelist.items()):
            input_data += convert_input_to_namelist_entry(key, value)
        input_data += '/\n'

        # Write CI cards now
        if manual_climbing_image:
            input_data += climbing_image_card

        if input_params:
            raise InputValidationError(
                'The following namelists are specified in input_params, but are not valid namelists for the current '
                f'type of calculation: {",".join(list(input_params.keys()))}'
            )

        return input_data

    def prepare_for_submission(self, folder):
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :class:`~aiida.common.datastructures.CalcInfo` instance.
        """
        # pylint: disable=too-many-branches,too-many-statements
        import numpy as np

        local_copy_list = []
        remote_copy_list = []
        remote_symlink_list = []

        # Convert settings dictionary to have uppercase keys, or create an empty one if none was given.
        if 'settings' in self.inputs:
            settings_dict = _uppercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
            if 'ADDITIONAL_RETRIEVE_LIST' in settings_dict:
                warnings.warn(
                    'The key `ADDITIONAL_RETRIEVE_LIST` in the settings input is deprecated and will be removed in '
                    'the future. Use the `CalcJob.metadata.options.additional_retrieve_list` input instead.',
                    AiidaDeprecationWarning
                )
        else:
            settings_dict = {}

        first_structure = self.inputs.first_structure
        last_structure = self.inputs.last_structure

        # Check that the first and last image have the same cell
        if abs(np.array(first_structure.cell) - np.array(last_structure.cell)).max() > 1.e-4:
            raise InputValidationError('Different cell in the fist and last image')

        # Check that the first and last image have the same number of sites
        if len(first_structure.sites) != len(last_structure.sites):
            raise InputValidationError('Different number of sites in the fist and last image')

        # Check that sites in the initial and final structure have the same kinds
        if first_structure.get_site_kindnames() != last_structure.get_site_kindnames():
            raise InputValidationError(
                'Mismatch between the kind names and/or order between '
                'the first and final image'
            )

        # Check that a pseudo potential was specified for each kind present in the `StructureData`
        # self.inputs.pw.pseudos is a plumpy.utils.AttributesFrozendict
        kindnames = [kind.name for kind in first_structure.kinds]
        if set(kindnames) != set(self.inputs.pw.pseudos.keys()):
            formatted_pseudos = ', '.join(list(self.inputs.pw.pseudos.keys()))
            formatted_kinds = ', '.join(list(kindnames))
            raise InputValidationError(
                'Mismatch between the defined pseudos and the list of kinds of the structure.\n'
                f'Pseudos: {formatted_pseudos};\nKinds: {formatted_kinds}'
            )

        ##############################
        # END OF INITIAL INPUT CHECK #
        ##############################

        # Create the subfolder that will contain the pseudopotentials
        folder.get_subfolder(self._PSEUDO_SUBFOLDER, create=True)
        # Create the subfolder for the output data (sometimes Quantum ESPRESSO codes crash if the folder does not exist)
        folder.get_subfolder(self._OUTPUT_SUBFOLDER, create=True)

        # We first prepare the NEB-specific input file.
        neb_input_filecontent = self._generate_input_files(self.inputs.parameters, settings_dict)
        with folder.open(self.inputs.metadata.options.input_filename, 'w') as handle:
            handle.write(neb_input_filecontent)

        # We now generate the PW input files for each input structure
        local_copy_pseudo_list = []
        for i, structure in enumerate([first_structure, last_structure]):
            # We need to a pass a copy of the settings_dict for each structure
            this_settings_dict = copy.deepcopy(settings_dict)
            pw_input_filecontent, this_local_copy_pseudo_list = PwCalculation._generate_PWCPinputdata(  # pylint: disable=protected-access
                self.inputs.pw.parameters, this_settings_dict, self.inputs.pw.pseudos, structure, self.inputs.pw.kpoints
            )
            local_copy_pseudo_list += this_local_copy_pseudo_list
            with folder.open(f'pw_{i + 1}.in', 'w') as handle:
                handle.write(pw_input_filecontent)

        # We need to pop the settings that were used in the PW calculations
        for key in list(settings_dict.keys()):
            if key not in list(this_settings_dict.keys()):
                settings_dict.pop(key)

        # We avoid to copy twice the same pseudopotential to the same filename
        local_copy_pseudo_list = set(local_copy_pseudo_list)
        # We check that two different pseudopotentials are not copied
        # with the same name (otherwise the first is overwritten)
        if len({filename for (uuid, filename, local_path) in local_copy_pseudo_list}) < len(local_copy_pseudo_list):
            raise InputValidationError('Same filename for two different pseudopotentials')

        local_copy_list += local_copy_pseudo_list

        # If present, add also the Van der Waals table to the pseudo dir. Note that the name of the table is not checked
        # but should be the one expected by Quantum ESPRESSO.
        vdw_table = self.inputs.get('pw.vdw_table', None)
        if vdw_table:
            local_copy_list.append(
                (vdw_table.uuid, vdw_table.filename, os.path.join(self._PSEUDO_SUBFOLDER, vdw_table.filename))
            )

        # operations for restart
        parent_calc_folder = self.inputs.get('parent_folder', None)
        symlink = settings_dict.pop('PARENT_FOLDER_SYMLINK', self._default_symlink_usage)  # a boolean
        if symlink:
            if parent_calc_folder is not None:
                # I put the symlink to the old parent ./out folder
                remote_symlink_list.append((
                    parent_calc_folder.computer.uuid,
                    os.path.join(parent_calc_folder.get_remote_path(), self._OUTPUT_SUBFOLDER,
                                 '*'),  # asterisk: make individual symlinks for each file
                    self._OUTPUT_SUBFOLDER
                ))
                # and to the old parent prefix.path
                remote_symlink_list.append((
                    parent_calc_folder.computer.uuid,
                    os.path.join(parent_calc_folder.get_remote_path(), f'{self._PREFIX}.path'), f'{self._PREFIX}.path'
                ))
        else:
            # copy remote output dir and .path file, if specified
            if parent_calc_folder is not None:
                remote_copy_list.append((
                    parent_calc_folder.computer.uuid,
                    os.path.join(parent_calc_folder.get_remote_path(), self._OUTPUT_SUBFOLDER,
                                 '*'), self._OUTPUT_SUBFOLDER
                ))
                # and copy the old parent prefix.path
                remote_copy_list.append((
                    parent_calc_folder.computer.uuid,
                    os.path.join(parent_calc_folder.get_remote_path(), f'{self._PREFIX}.path'), f'{self._PREFIX}.path'
                ))

        # here we may create an aiida.EXIT file
        create_exit_file = settings_dict.pop('ONLY_INITIALIZATION', False)
        if create_exit_file:
            exit_filename = f'{self._PREFIX}.EXIT'
            with folder.open(exit_filename, 'w') as handle:
                handle.write('\n')

        calcinfo = CalcInfo()
        codeinfo = CodeInfo()

        calcinfo.uuid = self.uuid
        cmdline_params = settings_dict.pop('CMDLINE', [])
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.remote_symlink_list = remote_symlink_list
        # In neb calculations there is no input read from standard input!!
        codeinfo.cmdline_params = (['-input_images', '2'] + list(cmdline_params))
        codeinfo.stdout_name = self.inputs.metadata.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid
        calcinfo.codes_info = [codeinfo]

        # Retrieve the output files and the xml files
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(self.inputs.metadata.options.output_filename)
        calcinfo.retrieve_list.append((
            os.path.join(self._OUTPUT_SUBFOLDER, self._PREFIX + '_*[0-9]', 'PW.out'),  # source relative path (globbing)
            '.',  # destination relative path
            2  # depth to preserve
        ))

        for xml_filepath in self.xml_filepaths:  # pylint: disable=not-an-iterable
            calcinfo.retrieve_list.append([xml_filepath, '.', 3])

        calcinfo.retrieve_list += settings_dict.pop('ADDITIONAL_RETRIEVE_LIST', [])
        calcinfo.retrieve_list += self._internal_retrieve_list

        # We might still have parser options in the settings dictionary: pop them.
        _pop_parser_options(self, settings_dict)

        if settings_dict:
            unknown_keys = ', '.join(list(settings_dict.keys()))
            raise InputValidationError(f'`settings` contained unexpected keys: {unknown_keys}')

        return calcinfo
