# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the pp.x code of Quantum ESPRESSO."""
import os
import warnings

from aiida import orm
from aiida.common import datastructures, exceptions
from aiida.common.warnings import AiidaDeprecationWarning

from aiida_quantumespresso.calculations import _lowercase_dict, _uppercase_dict
from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry

from .base import CalcJob


def validate_parameters(value, ctx=None):  # pylint: disable=unused-argument
    """Validate 'parameters' dict."""
    parameters = value.get_dict()

    try:
        plot_num = parameters['INPUTPP']['plot_num']
    except KeyError:
        return 'parameter `INPUTPP.plot_num` must be explicitly set'

    # Check that a valid plot type is requested
    if plot_num in range(23) and plot_num not in [14, 15, 16]:  # Must be integer in range 0-22, but not 14-16:
        parameters['INPUTPP']['plot_num'] = int(plot_num)  # If this test passes, we can safely cast to int
    else:
        return '`INTPUTPP.plot_num` must be an integer in the range 0-23 excluding [14, 15, 16]'

    try:
        dimension = parameters['PLOT']['iflag']
    except KeyError:
        return 'parameter `PLOT.iflag` must be explicitly set'

    # Check for valid plot dimension:
    if dimension in range(5):  # Must be in range 0-4:
        parameters['PLOT']['iflag'] = int(dimension)
    else:
        return '`PLOT.iflag` must be an integer in the range 0-4'


class PpCalculation(CalcJob):
    """`CalcJob` implementation for the pp.x code of Quantum ESPRESSO."""

    # Default name of the subfolder inside 'parent_folder' from which you want to copy the files, in case the
    # parent_folder is of type FolderData
    _INPUT_SUBFOLDER = './out/'

    # In the PW Calculation plugin, these folder names and prefixes are fixed, so import them to reduce maintenance
    # pylint: disable=protected-access
    from aiida_quantumespresso.calculations import BasePwCpInputGenerator
    _OUTPUT_SUBFOLDER = BasePwCpInputGenerator._OUTPUT_SUBFOLDER
    _PREFIX = BasePwCpInputGenerator._PREFIX
    _PSEUDO_SUBFOLDER = BasePwCpInputGenerator._PSEUDO_SUBFOLDER
    _DEFAULT_INPUT_FILE = BasePwCpInputGenerator._DEFAULT_INPUT_FILE
    _DEFAULT_OUTPUT_FILE = BasePwCpInputGenerator._DEFAULT_OUTPUT_FILE
    # pylint: enable=protected-access

    # Grid data output file from first stage of pp calculation
    _FILPLOT = 'aiida.filplot'
    # Grid data output in desired format
    _FILEOUT = 'aiida.fileout'

    _default_namelists = ['INPUTPP', 'PLOT']

    # Keywords that cannot be set by the user but will be set by the plugin
    _blocked_keywords = [('INPUTPP', 'outdir', _OUTPUT_SUBFOLDER), ('INPUTPP', 'prefix', _PREFIX),
                         ('INPUTPP', 'filplot', _FILPLOT), ('PLOT', 'fileout', _FILEOUT),
                         ('PLOT', 'output_format', None)]

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('parent_folder', valid_type=(orm.RemoteData, orm.FolderData), required=True,
            help='Output folder of a completed `PwCalculation`')
        spec.input('parameters', valid_type=orm.Dict, required=True, validator=validate_parameters,
            help='Use a node that specifies the input parameters for the namelists')
        spec.input('settings', valid_type=orm.Dict, required=False,
            help='Optional parameters to affect the way the calculation job is performed.')
        spec.input('metadata.options.input_filename', valid_type=str, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=str, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.parser_name', valid_type=str, default='quantumespresso.pp')
        spec.input('metadata.options.withmpi', valid_type=bool, default=True)
        spec.input('metadata.options.keep_plot_file', valid_type=bool, required=False)
        spec.input('metadata.options.keep_data_files', valid_type=bool, default=False)
        spec.input('metadata.options.parse_data_files', valid_type=bool, default=True)

        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('output_data', valid_type=orm.ArrayData)
        spec.output_namespace('output_data_multiple', valid_type=orm.ArrayData, dynamic=True)
        spec.default_output_node = 'output_parameters'

        # Standard exceptions
        spec.exit_code(301, 'ERROR_NO_RETRIEVED_TEMPORARY_FOLDER',
            message='The retrieved temporary folder could not be accessed.')
        spec.exit_code(302, 'ERROR_OUTPUT_STDOUT_MISSING',
            message='The retrieved folder did not contain the required stdout output file.')
        spec.exit_code(303, 'ERROR_OUTPUT_XML_MISSING',
            message='The parent folder did not contain the required XML output file.')
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(311, 'ERROR_OUTPUT_STDOUT_PARSE',
            message='The stdout output file could not be parsed.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete.')
        spec.exit_code(340, 'ERROR_OUT_OF_WALLTIME_INTERRUPTED',
            message='The calculation stopped prematurely because it ran out of walltime but the job was killed by the '
                    'scheduler before the files were safely written to disk for a potential restart.')
        spec.exit_code(350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION',
            message='The parser raised an unexpected exception: {exception}')

        # Output datafile related exceptions
        spec.exit_code(330, 'ERROR_OUTPUT_DATAFILE_MISSING',
            message='The formatted data output file `{filename}` was not present in the retrieved (temporary) folder.')
        spec.exit_code(331, 'ERROR_OUTPUT_DATAFILE_READ',
            message='The formatted data output file `{filename}` could not be read.')
        spec.exit_code(332, 'ERROR_UNSUPPORTED_DATAFILE_FORMAT',
            message='The data file format is not supported by the parser')
        spec.exit_code(333, 'ERROR_OUTPUT_DATAFILE_PARSE',
            message='The formatted data output file `{filename}` could not be parsed: {exception}')
        # yapf: enable

    def prepare_for_submission(self, folder):  # pylint: disable=too-many-branches,too-many-statements
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :class:`~aiida.common.datastructures.CalcInfo` instance.
        """

        # Put the first-level keys as uppercase (i.e., namelist and card names) and the second-level keys as lowercase
        parameters = _uppercase_dict(self.inputs.parameters.get_dict(), dict_name='parameters')
        parameters = {k: _lowercase_dict(v, dict_name=k) for k, v in parameters.items()}

        # Same for settings.
        if 'settings' in self.inputs:
            settings = _uppercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
        else:
            settings = {}

        # Set default values. NOTE: this is different from PW/CP
        for blocked in self._blocked_keywords:
            namelist = blocked[0].upper()
            key = blocked[1].lower()
            value = blocked[2]

            if namelist in parameters:
                if key in parameters[namelist]:
                    raise exceptions.InputValidationError(
                        f"You cannot specify explicitly the '{key}' key in the '{namelist}' namelist."
                    )
            else:
                parameters[namelist] = {}
            parameters[namelist][key] = value

        # Restrict the plot output to the file types that we want to be able to parse
        dimension_to_output_format = {
            0: 0,  # Spherical integration -> Gnuplot, 1D
            1: 0,  # 1D -> Gnuplot, 1D
            2: 7,  # 2D -> Gnuplot, 2D
            3: 6,  # 3D -> Gaussian cube
            4: 0,  # Polar on a sphere -> # Gnuplot, 1D
        }
        parameters['PLOT']['output_format'] = dimension_to_output_format[parameters['PLOT']['iflag']]

        namelists_toprint = self._default_namelists

        input_filename = self.inputs.metadata.options.input_filename
        with folder.open(input_filename, 'w') as infile:
            for namelist_name in namelists_toprint:
                infile.write(f'&{namelist_name}\n')
                # namelist content; set to {} if not present, so that we leave an empty namelist
                namelist = parameters.pop(namelist_name, {})
                for key, value in sorted(namelist.items()):
                    infile.write(convert_input_to_namelist_entry(key, value))
                infile.write('/\n')

        # Check for specified namelists that are not expected
        if parameters:
            raise exceptions.InputValidationError(
                'The following namelists are specified in parameters, but are not valid namelists for the current type '
                f'of calculation: {",".join(list(parameters.keys()))}'
            )

        remote_copy_list = []
        local_copy_list = []

        source = self.inputs.get('parent_folder', None)

        if isinstance(source, orm.RemoteData):
            dirpath = os.path.join(source.get_remote_path(), self._INPUT_SUBFOLDER)
            remote_copy_list.append((source.computer.uuid, dirpath, self._OUTPUT_SUBFOLDER))
            dirpath = os.path.join(source.get_remote_path(), self._PSEUDO_SUBFOLDER)
            remote_copy_list.append((source.computer.uuid, dirpath, self._PSEUDO_SUBFOLDER))
        elif isinstance(source, orm.FolderData):
            local_copy_list.append((source.uuid, self._OUTPUT_SUBFOLDER, self._OUTPUT_SUBFOLDER))
            local_copy_list.append((source.uuid, self._PSEUDO_SUBFOLDER, self._PSEUDO_SUBFOLDER))

        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = settings.pop('CMDLINE', [])
        codeinfo.stdin_name = self.inputs.metadata.options.input_filename
        codeinfo.stdout_name = self.inputs.metadata.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid

        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list

        # Retrieve by default the output file
        calcinfo.retrieve_list = [self.inputs.metadata.options.output_filename]
        calcinfo.retrieve_temporary_list = []

        # Depending on the `plot_num` and the corresponding parameters, more than one pair of `filplot` + `fileout`
        # files may be written. In that case, the data files will have `filplot` as a prefix with some suffix to
        # distinguish them from one another. The `fileout` filename will be the full data filename with the `fileout`
        # value as a suffix.
        retrieve_tuples = [self._FILEOUT, (f'{self._FILPLOT}*{self._FILEOUT}', '.', 0)]
        if 'keep_plot_file' in self.inputs.metadata.options:
            self.inputs.metadata.options.keep_data_files = self.inputs.metadata.options.keep_plot_file
            warnings.warn(
                "The input parameter 'keep_plot_file' is deprecated and will be removed in version 5.0.0. "
                "Please use 'keep_data_files' instead.", AiidaDeprecationWarning
            )
        if self.inputs.metadata.options.keep_data_files:
            calcinfo.retrieve_list.extend(retrieve_tuples)
        # If we do not want to parse the retrieved files, temporary retrieval is meaningless
        elif self.inputs.metadata.options.parse_data_files:
            calcinfo.retrieve_temporary_list.extend(retrieve_tuples)

        return calcinfo
