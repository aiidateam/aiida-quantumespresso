# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the pp.x code of Quantum ESPRESSO."""
from __future__ import absolute_import

import os
import six
from six.moves import range

from aiida import orm
from aiida.common import datastructures, exceptions

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
        # yapf: disable
        super(PpCalculation, cls).define(spec)
        spec.input('parent_folder', valid_type=(orm.RemoteData, orm.FolderData), required=True,
            help='Output folder of a completed `PwCalculation`')
        spec.input('parameters', valid_type=orm.Dict, required=True, validator=validate_parameters,
            help='Use a node that specifies the input parameters for the namelists')
        spec.input('metadata.options.input_filename', valid_type=six.string_types, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=six.string_types, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.parser_name', valid_type=six.string_types, default='quantumespresso.pp')
        spec.input('metadata.options.withmpi', valid_type=bool, default=True)
        spec.input('metadata.options.keep_plot_file', valid_type=bool, default=False)

        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('output_data', valid_type=orm.ArrayData)
        spec.default_output_node = 'output_parameters'

        # Standard exceptions
        spec.exit_code(300, 'ERROR_NO_RETRIEVED_FOLDER',
            message='The retrieved folder data node could not be accessed.')
        spec.exit_code(301, 'ERROR_NO_RETRIEVED_TEMPORARY_FOLDER',
            message='The retrieved temporary folder could not be accessed.')
        spec.exit_code(302, 'ERROR_OUTPUT_STDOUT_MISSING',
            message='The retrieved folder did not contain the required stdout output file.')
        spec.exit_code(303, 'ERROR_OUTPUT_XML_MISSING',
            message='The parent folder did not contain the required XML output file.')
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete.')
        spec.exit_code(340, 'ERROR_OUT_OF_WALLTIME_INTERRUPTED',
            message='The calculation stopped prematurely because it ran out of walltime but the job was killed by the '
                    'scheduler before the files were safely written to disk for a potential restart.')
        spec.exit_code(350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION',
            message='The parser raised an unexpected exception.')

        # Output datafile related exceptions
        spec.exit_code(330, 'ERROR_OUTPUT_DATAFILE_MISSING',
            message='The retrieved folder did not contain the required formatted data output file.')
        spec.exit_code(331, 'ERROR_OUTPUT_DATAFILE_READ',
            message='The formatted data output file could not be read.')
        spec.exit_code(332, 'ERROR_UNSUPPORTED_DATAFILE_FORMAT',
            message='The data file format is not supported by the parser')

    def prepare_for_submission(self, folder):  # pylint: disable=too-many-branches,too-many-statements
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """

        # Put the first-level keys as uppercase (i.e., namelist and card names) and the second-level keys as lowercase
        parameters = _uppercase_dict(self.inputs.parameters.get_dict(), dict_name='parameters')
        parameters = {k: _lowercase_dict(v, dict_name=k) for k, v in six.iteritems(parameters)}

        # Set default values. NOTE: this is different from PW/CP
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
                infile.write(u'&{0}\n'.format(namelist_name))
                # namelist content; set to {} if not present, so that we leave an empty namelist
                namelist = parameters.pop(namelist_name, {})
                for key, value in sorted(six.iteritems(namelist)):
                    infile.write(convert_input_to_namelist_entry(key, value))
                infile.write(u'/\n')

        # Check for specified namelists that are not expected
        if parameters:
            raise exceptions.InputValidationError(
                'The following namelists are specified in parameters, but are '
                'not valid namelists for the current type of calculation: '
                '{}'.format(','.join(list(parameters.keys()))))

        remote_copy_list = []
        local_copy_list = []

        # Copy remote output dir
        parent_calc_folder = self.inputs.get('parent_folder', None)
        if isinstance(parent_calc_folder, orm.RemoteData):
            remote_copy_list.append((
                parent_calc_folder.computer.uuid,
                os.path.join(parent_calc_folder.get_remote_path(), self._INPUT_SUBFOLDER),
                self._OUTPUT_SUBFOLDER
            ))
            remote_copy_list.append((
                parent_calc_folder.computer.uuid,
                os.path.join(parent_calc_folder.get_remote_path(), self._PSEUDO_SUBFOLDER),
                self._PSEUDO_SUBFOLDER
            ))
        elif isinstance(parent_calc_folder, orm.FolderData):
            for filename in parent_calc_folder.list_object_names():
                local_copy_list.append((
                    parent_calc_folder.uuid,
                    filename,
                    os.path.join(self._OUTPUT_SUBFOLDER, filename)
                ))
                local_copy_list.append((
                    parent_calc_folder.uuid,
                    filename,
                    os.path.join(self._PSEUDO_SUBFOLDER, filename)
                ))

        codeinfo = datastructures.CodeInfo()
        codeinfo.stdin_name = self.inputs.metadata.options.input_filename
        codeinfo.stdout_name = self.inputs.metadata.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid

        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list

        # Retrieve by default the output file and plot file
        calcinfo.retrieve_list = []
        calcinfo.retrieve_temporary_list = []
        calcinfo.retrieve_list.append(self.inputs.metadata.options.output_filename)
        if self.inputs.metadata.options.keep_plot_file:
            calcinfo.retrieve_list.append(self._FILEOUT)
        else:
            calcinfo.retrieve_temporary_list.append(self._FILEOUT)

        return calcinfo
