# -*- coding: utf-8 -*-
"""Plugin to create a Quantum Espresso epw.x input file."""
from __future__ import absolute_import

import os
import numpy
import six

from aiida import orm
from aiida.common import datastructures, exceptions
from aiida.engine import CalcJob

from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.calculations import _lowercase_dict, _uppercase_dict
from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry


class EpwCalculation(CalcJob):
    """`CalcJob` implementation for the epw.x code of Quantum ESPRESSO."""

    # Keywords that cannot be set by the user but will be set by the plugin
    _blocked_keywords = [('INPUTEPW', 'outdir'), ('INPUTEPW', 'iverbosity'), ('INPUTEPW', 'prefix'), 
                         ('INPUTEPW', 'dvscf_dir'),
                         ('INPUTEPW', 'amass'), ('INPUTEPW', 'nq1'), ('INPUTEPW', 'nq2'), ('INPUTEPW', 'nq3'),
                         ('INPUTEPW', 'nk1'), ('INPUTEPW', 'nk2'), ('INPUTEPW', 'nk3')]

    _use_kpoints = True

    _compulsory_namelists = ['INPUTEPW']

    # Default input and output files
    _PREFIX = 'aiida'
    _DEFAULT_INPUT_FILE = 'aiida.in'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'
    _OUTPUT_XML_TENSOR_FILE_NAME = 'tensors.xml'
    _OUTPUT_SUBFOLDER = './out/'
    _SAVE_PREFIX = '/save/'

    # Not using symlink in pw to allow multiple nscf to run on top of the same scf
    _default_symlink_usage = False

    @classmethod
    def define(cls, spec):
        # yapf: disable
        super(EpwCalculation, cls).define(spec)
        spec.input('metadata.options.input_filename', valid_type=six.string_types, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=six.string_types, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.parser_name', valid_type=six.string_types, default='quantumespresso.epw')
        spec.input('metadata.options.withmpi', valid_type=bool, default=True)
        spec.input('kpoints', valid_type=orm.KpointsData, help='kpoint mesh')
        spec.input('qpoints', valid_type=orm.KpointsData, help='qpoint mesh')
        spec.input('parameters', valid_type=orm.Dict, help='')
        spec.input('settings', valid_type=orm.Dict, required=False, help='')
        spec.input('parent_folder_nscf', valid_type=orm.RemoteData, help='the folder of a completed nscf `PwCalculation`')
        spec.input('parent_folder_ph', valid_type=orm.RemoteData, help='the folder of a completed `PhCalculation`')
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.default_output_node = 'output_parameters'

        # Unrecoverable errors: resources like the retrieved folder or its expected contents are missing
        spec.exit_code(200, 'ERROR_NO_RETRIEVED_FOLDER',
            message='The retrieved folder data node could not be accessed.')
        spec.exit_code(210, 'ERROR_OUTPUT_STDOUT_MISSING',
            message='The retrieved folder did not contain the required stdout output file.')

        # Unrecoverable errors: required retrieved files could not be read, parsed or are otherwise incomplete
        spec.exit_code(300, 'ERROR_OUTPUT_FILES',
            message='Both the stdout and XML output files could not be read or parsed.')
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(311, 'ERROR_OUTPUT_STDOUT_PARSE',
            message='The stdout output file could not be parsed.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete.')
        spec.exit_code(350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION',
            message='The parser raised an unexpected exception.')

        # Significant errors but calculation can be used to restart
        spec.exit_code(400, 'ERROR_OUT_OF_WALLTIME',
            message='The calculation stopped prematurely because it ran out of walltime.')
        spec.exit_code(410, 'ERROR_CONVERGENCE_NOT_REACHED',
            message='The minimization cycle did not reach self-consistency.')

    def prepare_for_submission(self, folder):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        # pylint: disable=too-many-statements,too-many-branches
        local_copy_list = []
        remote_copy_list = []
        remote_symlink_list = []

        if 'settings' in self.inputs:
            settings = _uppercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
        else:
            settings = {}

        # Copy nscf folder
        parent_folder_nscf = self.inputs.parent_folder_nscf
        parent_calcs_nscf = parent_folder_nscf.get_incoming(node_class=orm.CalcJobNode).all()

        if not parent_calcs_nscf:
            raise exceptions.NotExistent('parent_folder<{}> has no parent calculation'.format(parent_folder_nscf.pk))
        elif len(parent_calcs_nscf) > 1:
            raise exceptions.UniquenessError(
                'parent_folder<{}> has multiple parent calculations'.format(parent_folder_nscf.pk))

        parent_calc_nscf = parent_calcs_nscf[0].node

        # Also, the parent calculation must be on the same computer
        if not self.node.computer.uuid == parent_calc_nscf.computer.uuid:
            raise exceptions.InputValidationError(
                'Calculation has to be launched on the same computer as that of the parent: {}'.format(
                    parent_calc_nscf.computer.get_name()))

        # put by default, default_parent_output_folder = ./out
        try:
            default_parent_output_folder_nscf = parent_calc_nscf.process_class._OUTPUT_SUBFOLDER  # pylint: disable=protected-access
        except AttributeError:
            try:
                default_parent_output_folder_nscf = parent_calc_nscf._get_output_folder()  # pylint: disable=protected-access
            except AttributeError:
                raise exceptions.InputValidationError('parent calculation does not have a default output subfolder')
        parent_calc_out_subfolder_nscf = settings.pop('PARENT_CALC_OUT_SUBFOLDER_NSCF', default_parent_output_folder_nscf)

        # Now phonon folder 
        parent_folder_ph = self.inputs.parent_folder_ph
        parent_calcs_ph = parent_folder_ph.get_incoming(node_class=orm.CalcJobNode).all()

        if not parent_calcs_ph:
            raise exceptions.NotExistent('parent_folder<{}> has no parent calculation'.format(parent_folder_ph.pk))
        elif len(parent_calcs_ph) > 1:
            raise exceptions.UniquenessError(
                'parent_folder<{}> has multiple parent calculations'.format(parent_folder_ph.pk))

        parent_calc_ph = parent_calcs_ph[0].node

        # Also, the parent calculation must be on the same computer
        if not self.node.computer.uuid == parent_calc_ph.computer.uuid:
            raise exceptions.InputValidationError(
                'Calculation has to be launched on the same computer as that of the parent: {}'.format(
                    parent_calc_ph.computer.get_name()))

        # put by default, default_parent_output_folder = ./out
        try:
            default_parent_output_folder_ph = parent_calc_ph.process_class._OUTPUT_SUBFOLDER  # pylint: disable=protected-access
        except AttributeError:
            try:
                default_parent_output_folder_ph = parent_calc_ph._get_output_folder()  # pylint: disable=protected-access
            except AttributeError:
                raise exceptions.InputValidationError('parent calculation does not have a default output subfolder')
        parent_calc_out_subfolder_ph = settings.pop('PARENT_CALC_OUT_SUBFOLDER_PH', default_parent_output_folder_ph)




        # I put the first-level keys as uppercase (i.e., namelist and card names) and the second-level keys as lowercase
        parameters = _uppercase_dict(self.inputs.parameters.get_dict(), dict_name='parameters')
        parameters = {k: _lowercase_dict(v, dict_name=k) for k, v in six.iteritems(parameters)}

        if 'INPUTEPW' not in parameters:
            raise exceptions.InputValidationError('required namelist INPUTEPW not specified')

        parameters['INPUTEPW']['outdir'] = self._OUTPUT_SUBFOLDER
        parameters['INPUTEPW']['iverbosity'] = 1
        parameters['INPUTEPW']['prefix'] = self._PREFIX
        #parameters['INPUTPH']['fildyn'] = self._SAVE_PREFIX

        try:
            mesh, offset = self.inputs.qpoints.get_kpoints_mesh()

            if any([i != 0. for i in offset]):
                raise NotImplementedError(
                    'Computation of electron-phonon on a mesh with non zero offset is not implemented, at the level of epw.x')

            parameters['INPUTEPW']['nq1'] = mesh[0]
            parameters['INPUTEPW']['nq2'] = mesh[1]
            parameters['INPUTEPW']['nq3'] = mesh[2]

            postpend_text = None
        except:
            raise exceptions.InputValidationError('Cannot get the coarse q-point')

       	try:
            mesh, offset = self.inputs.kpoints.get_kpoints_mesh()

            if any([i != 0. for i in offset]):
                raise NotImplementedError(
                    'Computation of electron-phonon on a mesh with non zero offset is not implemented, at the level of epw.x')

            parameters['INPUTEPW']['nk1'] = mesh[0]
            parameters['INPUTEPW']['nk2'] = mesh[1]
            parameters['INPUTEPW']['nk3'] = mesh[2]

            postpend_text = None
        except:
            raise exceptions.InputValidationError('Cannot get the coarse k-point')


        # customized namelists, otherwise not present in the distributed epw code
        try:
            namelists_toprint = settings.pop('NAMELISTS')
            if not isinstance(namelists_toprint, list):
                raise exceptions.InputValidationError(
                    "The 'NAMELISTS' value, if specified in the settings input "
                    'node, must be a list of strings')
        except KeyError:  # list of namelists not specified in the settings; do automatic detection
            namelists_toprint = self._compulsory_namelists


        with folder.open(self.metadata.options.input_filename, 'w') as infile:
            for namelist_name in namelists_toprint:
                infile.write(u'&{0}\n'.format(namelist_name))
                # namelist content; set to {} if not present, so that we leave an empty namelist
                namelist = parameters.pop(namelist_name, {})
                for key, value in sorted(six.iteritems(namelist)):
                    infile.write(convert_input_to_namelist_entry(key, value))
                infile.write(u'/\n')

            # add list of qpoints if required
            if postpend_text is not None:
                infile.write(postpend_text)

        if parameters:
            raise exceptions.InputValidationError(
                'The following namelists are specified in parameters, but are '
                'not valid namelists for the current type of calculation: '
                '{}'.format(','.join(list(parameters.keys()))))

        # copy the parent scratch
        symlink = settings.pop('PARENT_FOLDER_SYMLINK', self._default_symlink_usage)  # a boolean
        if symlink:
            # I create a symlink to each file/folder in the parent ./out
            folder.get_subfolder(self._OUTPUT_SUBFOLDER, create=True)

            remote_symlink_list.append((
                parent_folder_nscf.computer.uuid,
                os.path.join(parent_folder_nscf.get_remote_path(), parent_calc_out_subfolder_nscf, '*'),
                self._OUTPUT_SUBFOLDER
            ))

            # I also create a symlink for the ./save folder
            #remote_symlink_list.append((
            #    parent_folder_ph.computer.uuid,
            #    os.path.join(parent_folder_ph.get_remote_path(), self._get_save_folder_ph()),
            #    self._get_save_folder()
            #))
        else:
            # here I copy the whole folder ./out
            remote_copy_list.append((
                parent_folder_nscf.computer.uuid,
                os.path.join(parent_folder_nscf.get_remote_path(), parent_calc_out_subfolder_nscf),
                self._OUTPUT_SUBFOLDER
            ))
            # I also copy the ./save folder
            #remote_copy_list.append((
            #    parent_folder_ph.computer.uuid,
            #    os.path.join(parent_folder_ph.get_remote_path(), self._get_pseudo_folder()),
            #    self._get_save_folder()
            #))


        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = (list(settings.pop('CMDLINE', [])) + ['-in', self.metadata.options.input_filename])
        codeinfo.stdout_name = self.metadata.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid

        calcinfo = datastructures.CalcInfo()
        calcinfo.uuid = str(self.uuid)
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.remote_symlink_list = remote_symlink_list

        # Retrieve by default the output file and the xml file
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(self.metadata.options.output_filename)

        if settings:
            unknown_keys = ', '.join(list(settings.keys()))
            raise exceptions.InputValidationError('`settings` contained unexpected keys: {}'.format(unknown_keys))

        return calcinfo

#    @staticmethod
#    def _get_pseudo_folder():
#        """Get the calculation-specific pseudo folder (relative path).
#
#        Default given by PwCalculation._PSEUDO_SUBFOLDER
#        """
#        return PwCalculation._PSEUDO_SUBFOLDER  # pylint: disable=protected-access
