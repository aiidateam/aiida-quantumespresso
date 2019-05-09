# -*- coding: utf-8 -*-
"""
Plugin to create a Quantum Espresso input file for a generic post-processing
(or similar) code that only requires a few namelists (plus possibly some text
afterwards).
"""
from __future__ import absolute_import
import os
from aiida.common import InputValidationError
from aiida.common import CalcInfo, CodeInfo
from aiida.orm import Dict, Code
from aiida.orm import RemoteData, FolderData, SinglefileData
from aiida.engine import CalcJob
from aiida.common import CodeInfo
from aiida.plugins import ParserFactory
from aiida_quantumespresso.calculations import _lowercase_dict, _uppercase_dict
from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry
import six

# TODO: review usage of class variables and input names
# TODO: add exit codes

class NamelistsCalculation(CalcJob):
    """
    Generic plugin to manage simple post-processing tools of the
    Quantum ESPRESSO distribution (http://www.quantum-espresso.org/)
    that accept as input a Fortran-style namelist.
    """

    # Default name of the subfolder inside 'parent_folder'
    # from which you want to copy the files, in case
    # the parent_folder is of type FolderData
    _INPUT_SUBFOLDER = "./out/"
    # Default name of the subfolder inside 'parent_folder'
    # from which you want to copy the files, in case
    # the parent_folder is of type RemoteData,
    # unless the user specified a SETTINGS->parent_calc_out_subfolder
    # value
    _default_parent_output_folder = './out/'
    # Default name of the subfolder that you want to create
    # in the working directory, intended to contain the outputs
    # of the code, and in which you want to place the files
    # taken from parent_folder/INPUT_SUBFOLDER, in case
    # the parent_folder is of type RemoteData or FolderData
    _OUTPUT_SUBFOLDER = './out/'
    _PREFIX = 'aiida'
    #_INPUT_FILE_NAME = 'aiida.in'
    #_OUTPUT_FILE_NAME = 'aiida.out'
    _internal_retrieve_list = []
    _default_namelists = ['INPUTPP']
    _blocked_keywords = [] # a list of tuples with key and value fixed
    #_parent_folder_type = (RemoteData, FolderData, SinglefileData)
    _retrieve_singlefile_list = []
    # Default input and output files
    _DEFAULT_INPUT_FILE = 'aiida.in'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'
    # Default parser
    _default_parser = None

    @classmethod
    def define(cls, spec):
        super(NamelistsCalculation, cls).define(spec)
        spec.input('code', valid_type=Code, help='')
        spec.input('metadata.options.input_filename', valid_type=six.string_types, default=cls._DEFAULT_INPUT_FILE, non_db=True)
        spec.input('metadata.options.output_filename', valid_type=six.string_types, default=cls._DEFAULT_OUTPUT_FILE, non_db=True)
        spec.input('metadata.options.parser_name', valid_type=six.string_types, default=cls._default_parser, non_db=True)
        spec.input('metadata.options.withmpi', valid_type=bool, default=True)  # Override default value False
        spec.input('parameters', valid_type=Dict, required=False, help='Use a node that specifies the input parameters for the namelists')
        spec.input('settings', valid_type=Dict, required=False, help='Use an additional node for special settings')
        spec.input('parent_folder', valid_type=(RemoteData, FolderData, SinglefileData), required=False, help='Use a local or remote folder as parent folder (for restarts and similar)')

    def _get_following_text(self, settings_dict):
        """
        By default, no text follows the namelists section. In case any additional information
        needs to be added to the input file, a subclass can override this method.
        If this function needs any data from the settings input, it should pop() them from settings_dict,
        so that prepare_for_submission() can check if any keys are left and are therefore unrecognized.
        """
        return u''

    def prepare_for_submission(self, folder):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        local_copy_list = []
        remote_copy_list = []

        # Settings converted to uppercase
        if 'settings' in self.inputs:
            settings_dict = _uppercase_dict(self.inputs.settings.get_dict(),
                                            dict_name='settings')
        else:
            settings_dict = {}

        following_text = self._get_following_text(settings_dict)

        ##############################
        # END OF INITIAL INPUT CHECK #
        ##############################

        # I put the first-level keys as uppercase (i.e., namelist and card names)
        # and the second-level keys as lowercase
        # (deeper levels are unchanged)
        if 'parameters' in self.inputs:
            input_params = _uppercase_dict(self.inputs.parameters.get_dict(),
                                           dict_name='parameters')
            input_params = {k: _lowercase_dict(v, dict_name=k)
                            for k, v in six.iteritems(input_params)}
        else:
            input_params = {}

        # set default values. NOTE: this is different from PW/CP
        for blocked in self._blocked_keywords:
            namelist = blocked[0].upper()
            key = blocked[1].lower()
            value = blocked[2]

            if namelist in input_params:
                if key in input_params[namelist]:
                    raise InputValidationError(
                        "You cannot specify explicitly the '{}' key in the '{}' "
                        "namelist.".format(key, namelist))
            else:
                input_params[namelist] = {}
            input_params[namelist][key] = value

        # =================== NAMELISTS AND CARDS ========================
        try:
            namelists_toprint = settings_dict.pop('NAMELISTS')
            if not isinstance(namelists_toprint, list):
                raise InputValidationError(
                    "The 'NAMELISTS' value, if specified in the settings input "
                    "node, must be a list of strings")
        except KeyError: # list of namelists not specified; do automatic detection
            namelists_toprint = self._default_namelists

        input_filename = self.inputs.metadata.options.input_filename
        with folder.open(input_filename,'w') as infile:
            for namelist_name in namelists_toprint:
                infile.write(u"&{0}\n".format(namelist_name))
                # namelist content; set to {} if not present, so that we leave an
                # empty namelist
                namelist = input_params.pop(namelist_name, {})
                for k, v in sorted(six.iteritems(namelist)):
                    infile.write(convert_input_to_namelist_entry(k,v))
                infile.write(u"/\n")

            # Write remaning text now, if any
            infile.write(following_text)

        # Check for specified namelists that are not expected
        if input_params:
            raise InputValidationError(
                "The following namelists are specified in input_params, but are "
                "not valid namelists for the current type of calculation: "
                "{}".format(",".join(list(input_params.keys()))))

        # copy remote output dir, if specified
        parent_calc_folder =  self.inputs.get('parent_folder', None)
        if parent_calc_folder is not None:
            if isinstance(parent_calc_folder, RemoteData):
                parent_calc_out_subfolder = settings_dict.pop('PARENT_CALC_OUT_SUBFOLDER',
                                              self._INPUT_SUBFOLDER)
                remote_copy_list.append(
                         (parent_calc_folder.computer.uuid,
                          os.path.join(parent_calc_folder.get_remote_path(),
                                       parent_calc_out_subfolder),
                          self._OUTPUT_SUBFOLDER))
            elif isinstance(parent_calc_folder, FolderData):
                # TODO: test me, especially with deep relative paths.
                for filename in parent_calc_folder.list_object_names():
                    local_copy_list.append(
                        (parent_calc_folder.uuid,
                         filename,
                         os.path.join(self._OUTPUT_SUBFOLDER, filename))
                    )
            elif isinstance(parent_calc_folder, SinglefileData):
                # TODO: test me
                single_file =  parent_calc_folder
                local_copy_list.append(
                    (single_file.uuid, single_file.filename, single_file.filename)
                 )

        calcinfo = CalcInfo()

        calcinfo.uuid = self.uuid
        # Empty command line by default
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list

        codeinfo = CodeInfo()
        codeinfo.cmdline_params = settings_dict.pop('CMDLINE', [])
        codeinfo.stdin_name = self.inputs.metadata.options.input_filename
        codeinfo.stdout_name = self.inputs.metadata.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid
        calcinfo.codes_info = [codeinfo]

        # Retrieve by default the output file and the xml file
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(self.inputs.metadata.options.output_filename)
        settings_retrieve_list = settings_dict.pop('ADDITIONAL_RETRIEVE_LIST', [])
        calcinfo.retrieve_list += settings_retrieve_list
        calcinfo.retrieve_list += self._internal_retrieve_list

        calcinfo.retrieve_singlefile_list = self._retrieve_singlefile_list

        # We might still have parser options in the settings dictionary: pop them
        # We need an instance of the parser class to get the parser options key (typically 'parser_options')
        parser_name = self.inputs.get('metadata.options.parser_name', None)
        if parser_name is not None:
            Parserclass = ParserFactory(parser_name)
            parser = Parserclass(self)
            try:
                parser_opts = parser.get_parser_settings_key().upper()
                settings_dict.pop(parser_opts)
            except (KeyError, AttributeError):
                # the key parser_opts isn't inside the dictionary,
                # or the parser doesn't have a method get_parser_settings_key().
                pass

        if settings_dict:
            raise InputValidationError("The following keys have been found in "
                  "the settings input node, but were not understood: {}".format(
                  ",".join(list(settings_dict.keys()))))

        return calcinfo
