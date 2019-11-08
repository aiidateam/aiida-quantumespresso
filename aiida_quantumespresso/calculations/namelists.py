# -*- coding: utf-8 -*-
"""Plugin to create a Quantum Espresso input file for a generic post-processing (or similar) code that only requires a
few namelists (plus possibly some text afterwards)."""
from __future__ import absolute_import

import os
import six

from aiida.common import datastructures, exceptions
from aiida.engine import CalcJob
from aiida.orm import Dict
from aiida.orm import RemoteData, FolderData, SinglefileData

from aiida_quantumespresso.calculations import _lowercase_dict, _uppercase_dict, _pop_parser_options
from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry


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
        # yapf: disable
        super(NamelistsCalculation, cls).define(spec)
        spec.input('metadata.options.input_filename', valid_type=six.string_types, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=six.string_types, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.withmpi', valid_type=bool, default=True)  # Override default value False
        if cls._default_parser is not None:
            spec.input('metadata.options.parser_name', valid_type=six.string_types, default=cls._default_parser)
        spec.input('parameters', valid_type=Dict, required=False,
            help='Use a node that specifies the input parameters for the namelists')
        spec.input('settings', valid_type=Dict, required=False,
            help='Use an additional node for special settings')
        spec.input('parent_folder', valid_type=(RemoteData, FolderData, SinglefileData), required=False,
            help='Use a local or remote folder as parent folder (for restarts and similar)')

    def _get_following_text(self):
        """Return any optional text that is to be written after the normal namelists in the input file.

        By default, no text follows the namelists section. If in a sub class, any additional information needs to be
        added to the input file, this method can be overridden to return the lines that should be appended.
        """
        # pylint: disable=no-self-use
        return u''

    def prepare_for_submission(self, folder):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
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
            parameters = {k: _lowercase_dict(v, dict_name=k) for k, v in six.iteritems(parameters)}
        else:
            parameters = {}

        # Force default values for blocked keywords. NOTE: this is different from PW/CP
        for blocked in self._blocked_keywords:
            namelist = blocked[0].upper()
            key = blocked[1].lower()
            value = blocked[2]
            if namelist in parameters:
                if key in parameters[namelist]:
                    raise exceptions.InputValidationError(
                        "You cannot specify explicitly the '{}' key in the '{}' "
                        'namelist.'.format(key, namelist))
            else:
                parameters[namelist] = {}
            parameters[namelist][key] = value

        # =================== NAMELISTS AND CARDS ========================
        try:
            namelists_toprint = settings.pop('NAMELISTS')
            if not isinstance(namelists_toprint, list):
                raise exceptions.InputValidationError(
                    "The 'NAMELISTS' value, if specified in the settings input node, must be a list of strings")
        except KeyError:  # list of namelists not specified; do automatic detection
            namelists_toprint = self._default_namelists

        input_filename = self.inputs.metadata.options.input_filename
        with folder.open(input_filename, 'w') as infile:
            for namelist_name in namelists_toprint:
                infile.write(u'&{0}\n'.format(namelist_name))
                # namelist content; set to {} if not present, so that we leave an empty namelist
                namelist = parameters.pop(namelist_name, {})
                for key, value in sorted(six.iteritems(namelist)):
                    infile.write(convert_input_to_namelist_entry(key, value))
                infile.write(u'/\n')

            # Write remaning text now, if any
            infile.write(following_text)

        # Check for specified namelists that are not expected
        if parameters:
            raise exceptions.InputValidationError(
                'The following namelists are specified in parameters, but are '
                'not valid namelists for the current type of calculation: '
                '{}'.format(','.join(list(parameters.keys()))))

        remote_copy_list = []
        local_copy_list = []

        # copy remote output dir, if specified
        parent_calc_folder = self.inputs.get('parent_folder', None)
        if parent_calc_folder is not None:
            if isinstance(parent_calc_folder, RemoteData):
                parent_calc_out_subfolder = settings.pop('PARENT_CALC_OUT_SUBFOLDER', self._INPUT_SUBFOLDER)
                remote_copy_list.append((
                    parent_calc_folder.computer.uuid,
                    os.path.join(parent_calc_folder.get_remote_path(), parent_calc_out_subfolder),
                    self._OUTPUT_SUBFOLDER
                ))
            elif isinstance(parent_calc_folder, FolderData):
                for filename in parent_calc_folder.list_object_names():
                    local_copy_list.append((
                        parent_calc_folder.uuid,
                        filename,
                        os.path.join(self._OUTPUT_SUBFOLDER, filename)
                    ))
            elif isinstance(parent_calc_folder, SinglefileData):
                single_file = parent_calc_folder
                local_copy_list.append(
                    (single_file.uuid, single_file.filename, single_file.filename)
                )

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
            raise exceptions.InputValidationError('`settings` contained unexpected keys: {}'.format(unknown_keys))

        return calcinfo
