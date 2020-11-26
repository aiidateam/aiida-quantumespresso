# -*- coding: utf-8 -*-
"""Plugin to create a Quantum Espresso input file for a generic post-processing code.

These codes typically only require a few namelists (plus possibly some text afterwards).
"""
import os

from aiida.common import datastructures, exceptions
from aiida.orm import Dict
from aiida.orm import RemoteData, FolderData, SinglefileData

from aiida_quantumespresso.calculations import _lowercase_dict, _uppercase_dict, _pop_parser_options
from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry

from .base import CalcJob


class NamelistsCalculation(CalcJob):
    """`CalcJob` implementation to serve as base class for simple post-processing tools of Quantum ESPRESSO."""

    # Default name of the subfolder inside 'parent_folder' from which you want to copy the files, in case the
    # parent_folder is of type FolderData
    _INPUT_SUBFOLDER = './out/'

    # Default name of the subfolder inside 'parent_folder' from which you want to copy the files, in case the
    # parent_folder is of type RemoteData, unless the user specified `parent_calc_out_subfolder` in the `settings`
    _default_parent_output_folder = './out/'

    # Default name of the subfolder that you want to create in the working directory, intended to contain the outputs
    # of the code, and in which you want to place the files taken from parent_folder/INPUT_SUBFOLDER, in case the
    # parent_folder is of type RemoteData or FolderData
    _OUTPUT_SUBFOLDER = './out/'
    _PREFIX = 'aiida'
    _internal_retrieve_list = []
    _default_namelists = ['INPUTPP']
    _blocked_keywords = []  # a list of tuples with key and value fixed

    _retrieve_singlefile_list = []

    _DEFAULT_INPUT_FILE = 'aiida.in'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'

    _default_parser = None

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('metadata.options.input_filename', valid_type=str, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=str, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.withmpi', valid_type=bool, default=True)  # Override default value False
        if cls._default_parser is not None:
            spec.input('metadata.options.parser_name', valid_type=str, default=cls._default_parser)
        spec.input('parameters', valid_type=Dict, required=False,
            help='Use a node that specifies the input parameters for the namelists')
        spec.input('settings', valid_type=Dict, required=False,
            help='Use an additional node for special settings')
        spec.input('parent_folder', valid_type=(RemoteData, FolderData, SinglefileData), required=False,
            help='Use a local or remote folder as parent folder (for restarts and similar)')
        # yapf: enable

    def _get_following_text(self):
        """Return any optional text that is to be written after the normal namelists in the input file.

        By default, no text follows the namelists section. If in a sub class, any additional information needs to be
        added to the input file, this method can be overridden to return the lines that should be appended.
        """
        # pylint: disable=no-self-use
        return ''

    @classmethod
    def set_blocked_keywords(cls, parameters):
        """Force default values for blocked keywords. NOTE: this is different from PW/CP."""
        for blocked in cls._blocked_keywords:
            namelist = blocked[0].upper()
            key = blocked[1].lower()
            value = blocked[2]
            if namelist in parameters:
                if key in parameters[namelist]:
                    raise exceptions.InputValidationError(
                        "You cannot specify explicitly the '{}' key in the '{}' "
                        'namelist.'.format(key, namelist)
                    )
            else:
                parameters[namelist] = {}
            parameters[namelist][key] = value

        return parameters

    @staticmethod
    def filter_namelists(parameters, namelists_toprint, check_remaining=True):
        """Select only the given namelists from a parameter dict.

        :param parameters: 'dict' containing the fortran namelists and parameters to be used.
        :param namelists_toprint:  'list' containing the namelists to be selected from 'parameters'.
          If a given namelist is not present, an empty dict is used
        :param check_remaining: 'bool', if True, raise an exception if more namelists other than
          the ones given in 'namelist_toprint' are present in 'parameters'.
        :return: 'dict' of namelists.
        :raise: InputValidationError
        """
        filtered = {}
        for key in namelists_toprint:
            # namelist content; set to {} if not present, so that we leave an empty namelist
            filtered[key] = parameters.pop(key, {})

        if parameters and check_remaining:
            raise exceptions.InputValidationError(
                'The following namelists are specified in parameters, but are '
                'not valid namelists for the current type of calculation: '
                '{}'.format(','.join(list(parameters.keys())))
            )

        return filtered

    @staticmethod
    def generate_input_file(parameters):
        """Generate namelist input_file content given a dict of parameters.

        :param parameters: 'dict' containing the fortran namelists and parameters to be used.
          e.g.: {'CONTROL':{'calculation':'scf'}, 'SYSTEM':{'ecutwfc':30}}
        :return: 'str' containing the input_file content a plain text.
        """

        file_lines = []
        for namelist_name, namelist in parameters.items():
            file_lines.append(f'&{namelist_name}')
            for key, value in sorted(namelist.items()):
                file_lines.append(convert_input_to_namelist_entry(key, value)[:-1])
            file_lines.append('/')

        return '\n'.join(file_lines)

    def prepare_for_submission(self, folder):
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :py:`~aiida.common.datastructures.CalcInfo` instance.
        """
        # pylint: disable=too-many-branches,too-many-statements
        if 'settings' in self.inputs:
            settings = _uppercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
        else:
            settings = {}

        following_text = self._get_following_text()

        # Put the first-level keys as uppercase (i.e., namelist and card names) and the second-level keys as lowercase
        if 'parameters' in self.inputs:
            parameters = _uppercase_dict(self.inputs.parameters.get_dict(), dict_name='parameters')
            parameters = {k: _lowercase_dict(v, dict_name=k) for k, v in parameters.items()}
        else:
            parameters = {}

        # =================== NAMELISTS AND CARDS ========================
        try:
            namelists_toprint = settings.pop('NAMELISTS')
            if not isinstance(namelists_toprint, list):
                raise exceptions.InputValidationError(
                    "The 'NAMELISTS' value, if specified in the settings input node, must be a list of strings"
                )
        except KeyError:  # list of namelists not specified; do automatic detection
            namelists_toprint = self._default_namelists

        parameters = self.set_blocked_keywords(parameters)
        parameters = self.filter_namelists(parameters, namelists_toprint)
        file_content = self.generate_input_file(parameters)
        file_content += '\n' + following_text
        input_filename = self.inputs.metadata.options.input_filename
        with folder.open(input_filename, 'w') as infile:
            infile.write(file_content)

        symlink = settings.pop('PARENT_FOLDER_SYMLINK', False)

        remote_copy_list = []
        local_copy_list = []
        remote_symlink_list = []

        ptr = remote_symlink_list if symlink else remote_copy_list

        # copy remote output dir, if specified
        parent_calc_folder = self.inputs.get('parent_folder', None)
        if parent_calc_folder is not None:
            if isinstance(parent_calc_folder, RemoteData):
                parent_calc_out_subfolder = settings.pop('PARENT_CALC_OUT_SUBFOLDER', self._INPUT_SUBFOLDER)
                ptr.append((
                    parent_calc_folder.computer.uuid,
                    os.path.join(parent_calc_folder.get_remote_path(),
                                 parent_calc_out_subfolder), self._OUTPUT_SUBFOLDER
                ))
            elif isinstance(parent_calc_folder, FolderData):
                for filename in parent_calc_folder.list_object_names():
                    local_copy_list.append(
                        (parent_calc_folder.uuid, filename, os.path.join(self._OUTPUT_SUBFOLDER, filename))
                    )
            elif isinstance(parent_calc_folder, SinglefileData):
                single_file = parent_calc_folder
                local_copy_list.append((single_file.uuid, single_file.filename, single_file.filename))

        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = settings.pop('CMDLINE', [])
        codeinfo.stdin_name = self.inputs.metadata.options.input_filename
        codeinfo.stdout_name = self.inputs.metadata.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid

        calcinfo = datastructures.CalcInfo()
        calcinfo.uuid = str(self.uuid)
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.remote_symlink_list = remote_symlink_list

        # Retrieve by default the output file and the xml file
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(self.inputs.metadata.options.output_filename)
        settings_retrieve_list = settings.pop('ADDITIONAL_RETRIEVE_LIST', [])
        calcinfo.retrieve_list += settings_retrieve_list
        calcinfo.retrieve_list += self._internal_retrieve_list

        calcinfo.retrieve_singlefile_list = self._retrieve_singlefile_list

        # We might still have parser options in the settings dictionary: pop them.
        _pop_parser_options(self, settings)

        if settings:
            unknown_keys = ', '.join(list(settings.keys()))
            raise exceptions.InputValidationError(f'`settings` contained unexpected keys: {unknown_keys}')

        return calcinfo
