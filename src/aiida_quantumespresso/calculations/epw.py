# -*- coding: utf-8 -*-
"""Plugin to create a Quantum Espresso epw.x input file."""
import os

from aiida import orm
from aiida.common import datastructures, exceptions
import numpy as np

from aiida_quantumespresso.calculations import _lowercase_dict, _uppercase_dict
from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry

from .base import CalcJob


class EpwCalculation(CalcJob):
    """`CalcJob` implementation for the epw.x code of Quantum ESPRESSO."""

    # Keywords that cannot be set by the user but will be set by the plugin
    _blocked_keywords = [('INPUTEPW', 'outdir'), ('INPUTEPW', 'verbosity'), ('INPUTEPW', 'prefix'),
                         ('INPUTEPW', 'dvscf_dir'), ('INPUTEPW', 'amass'), ('INPUTEPW', 'nq1'), ('INPUTEPW', 'nq2'),
                         ('INPUTEPW', 'nq3'), ('INPUTEPW', 'nk1'), ('INPUTEPW', 'nk2'), ('INPUTEPW', 'nk3')]

    _use_kpoints = True

    _compulsory_namelists = ['INPUTEPW']

    # Default input and output files
    _PREFIX = 'aiida'
    _DEFAULT_INPUT_FILE = 'aiida.in'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'
    _OUTPUT_XML_TENSOR_FILE_NAME = 'tensors.xml'
    _OUTPUT_SUBFOLDER = './out/'
    _SAVE_PREFIX = '/save/'
    _FOLDER_SAVE = 'save'
    _VERBOSITY = 'high'
    _FOLDER_DYNAMICAL_MATRIX = 'DYN_MAT'

    # Not using symlink in pw to allow multiple nscf to run on top of the same scf
    _default_symlink_usage = False

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('metadata.options.input_filename', valid_type=str, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=str, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.withmpi', valid_type=bool, default=True)
        spec.input('kpoints', valid_type=orm.KpointsData, help='coarse kpoint mesh')
        spec.input('qpoints', valid_type=orm.KpointsData, help='coarse qpoint mesh')
        spec.input('kfpoints', valid_type=orm.KpointsData, help='fine kpoint mesh')
        spec.input('qfpoints', valid_type=orm.KpointsData, help='fine qpoint mesh')
        spec.input('parameters', valid_type=orm.Dict, help='')
        spec.input('settings', valid_type=orm.Dict, required=False, help='')
        spec.input('parent_folder_nscf', valid_type=orm.RemoteData,
                 help='the folder of a completed nscf `PwCalculation`')
        spec.input('parent_folder_ph', valid_type=orm.RemoteData, help='the folder of a completed `PhCalculation`')
        # yapf: enable

    def prepare_for_submission(self, folder):  # pylint: disable=too-many-statements,too-many-branches
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :class:`~aiida.common.datastructures.CalcInfo` instance.
        """

        def test_offset(offset):
            """Check if the grid has an offset."""
            if any(i != 0. for i in offset):
                raise NotImplementedError(
                    'Computation of electron-phonon on a mesh with non zero offset is not implemented, '
                    'at the level of epw.x'
                )

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
        parent_calc_nscf = parent_folder_nscf.creator

        if parent_calc_nscf is None:
            raise exceptions.NotExistent(f'parent_folder<{parent_folder_nscf.pk}> has no parent calculation')

        # Also, the parent calculation must be on the same computer
        if not self.node.computer.uuid == parent_calc_nscf.computer.uuid:
            computer_label = parent_calc_nscf.computer.get_name()
            raise exceptions.InputValidationError(
                f'Calculation has to be launched on the same computer as that of the parent: {computer_label}'
            )

        # put by default, default_parent_output_folder = ./out
        parent_calc_out_subfolder_nscf = parent_calc_nscf.process_class._OUTPUT_SUBFOLDER  # pylint: disable=protected-access

        # Now phonon folder
        parent_folder_ph = self.inputs.parent_folder_ph
        parent_calc_ph = parent_folder_ph.creator

        # Also, the parent calculation must be on the same computer
        if not self.node.computer.uuid == parent_calc_ph.computer.uuid:
            computer_label = parent_calc_nscf.computer.get_name()
            raise exceptions.InputValidationError(
                f'Calculation has to be launched on the same computer as that of the parent: {computer_label}'
            )

        # I put the first-level keys as uppercase (i.e., namelist and card names) and the second-level keys as lowercase
        parameters = _uppercase_dict(self.inputs.parameters.get_dict(), dict_name='parameters')
        parameters = {k: _lowercase_dict(v, dict_name=k) for k, v in parameters.items()}

        if 'INPUTEPW' not in parameters:
            raise exceptions.InputValidationError('required namelist INPUTEPW not specified')

        parameters['INPUTEPW']['outdir'] = self._OUTPUT_SUBFOLDER
        parameters['INPUTEPW']['verbosity'] = self._VERBOSITY
        parameters['INPUTEPW']['prefix'] = self._PREFIX

        try:
            mesh, offset = self.inputs.qpoints.get_kpoints_mesh()
            test_offset(offset)
            parameters['INPUTEPW']['nq1'] = mesh[0]
            parameters['INPUTEPW']['nq2'] = mesh[1]
            parameters['INPUTEPW']['nq3'] = mesh[2]
            postpend_text = None
        except NotImplementedError as exception:
            raise exceptions.InputValidationError('Cannot get the coarse q-point grid') from exception

        try:
            mesh, offset = self.inputs.kpoints.get_kpoints_mesh()
            test_offset(offset)
            parameters['INPUTEPW']['nk1'] = mesh[0]
            parameters['INPUTEPW']['nk2'] = mesh[1]
            parameters['INPUTEPW']['nk3'] = mesh[2]
            postpend_text = None
        except NotImplementedError as exception:
            raise exceptions.InputValidationError('Cannot get the coarse k-point grid') from exception

        try:
            mesh, offset = self.inputs.qfpoints.get_kpoints_mesh()
            test_offset(offset)
            parameters['INPUTEPW']['nqf1'] = mesh[0]
            parameters['INPUTEPW']['nqf2'] = mesh[1]
            parameters['INPUTEPW']['nqf3'] = mesh[2]
            postpend_text = None
        except NotImplementedError as exception:
            raise exceptions.InputValidationError('Cannot get the fine q-point grid') from exception

        try:
            mesh, offset = self.inputs.kfpoints.get_kpoints_mesh()
            test_offset(offset)
            parameters['INPUTEPW']['nkf1'] = mesh[0]
            parameters['INPUTEPW']['nkf2'] = mesh[1]
            parameters['INPUTEPW']['nkf3'] = mesh[2]
            postpend_text = None
        except NotImplementedError as exception:
            raise exceptions.InputValidationError('Cannot get the fine k-point grid') from exception

        # customized namelists, otherwise not present in the distributed epw code
        try:
            namelists_toprint = settings.pop('NAMELISTS')
            if not isinstance(namelists_toprint, list):
                raise exceptions.InputValidationError(
                    "The 'NAMELISTS' value, if specified in the settings input "
                    'node, must be a list of strings'
                )
        except KeyError:  # list of namelists not specified in the settings; do automatic detection
            namelists_toprint = self._compulsory_namelists

        # create the save folder with dvscf and dyn files.
        folder.get_subfolder(self._FOLDER_SAVE, create=True)

        # List of IBZ q-point to be added below EPW. To be removed when removed from EPW.
        qibz_ar = []
        for key, value in sorted(parent_folder_ph.creator.outputs.output_parameters.get_dict().items()):
            if key.startswith('dynamical_matrix_'):
                qibz_ar.append(value['q_point'])

        qibz_node = orm.ArrayData()
        qibz_node.set_array('qibz', np.array(qibz_ar))

        list_of_points = qibz_node.get_array('qibz')
        # Number of q-point in the irreducible Brillouin Zone.
        nqpt = len(list_of_points[0, :])

        # add here the list of point coordinates
        if len(list_of_points) > 1:
            postpend_text = f'{len(list_of_points)} cartesian\n'
            for points in list_of_points:
                postpend_text += '{0:18.10f} {1:18.10f} {2:18.10f} \n'.format(*points)  # pylint: disable=consider-using-f-string

        with folder.open(self.metadata.options.input_filename, 'w') as infile:
            for namelist_name in namelists_toprint:
                infile.write(f'&{namelist_name}\n')
                # namelist content; set to {} if not present, so that we leave an empty namelist
                namelist = parameters.pop(namelist_name, {})
                for key, value in sorted(namelist.items()):
                    infile.write(convert_input_to_namelist_entry(key, value))
                infile.write('/\n')

            # add list of qpoints if required
            if postpend_text is not None:
                infile.write(postpend_text)

        if parameters:
            raise exceptions.InputValidationError(
                'The following namelists are specified in parameters, but are not valid namelists for the current type '
                f'of calculation: {",".join(list(parameters.keys()))}'
            )

        # copy the parent scratch
        symlink = settings.pop('PARENT_FOLDER_SYMLINK', self._default_symlink_usage)  # a boolean
        if symlink:
            # I create a symlink to each file/folder in the parent ./out
            folder.get_subfolder(self._OUTPUT_SUBFOLDER, create=True)

            remote_symlink_list.append((
                parent_folder_nscf.computer.uuid,
                os.path.join(parent_folder_nscf.get_remote_path(), parent_calc_out_subfolder_nscf,
                             '*'), self._OUTPUT_SUBFOLDER
            ))

        else:
            # here I copy the whole folder ./out
            remote_copy_list.append((
                parent_folder_nscf.computer.uuid,
                os.path.join(parent_folder_nscf.get_remote_path(),
                             parent_calc_out_subfolder_nscf), self._OUTPUT_SUBFOLDER
            ))

        prefix = self._PREFIX

        for iqpt in range(1, nqpt + 1):
            label = str(iqpt)
            tmp_path = os.path.join(self._FOLDER_DYNAMICAL_MATRIX, 'dynamical-matrix-0')
            remote_copy_list.append((
                parent_folder_ph.computer.uuid, os.path.join(parent_folder_ph.get_remote_path(),
                                                             tmp_path), 'save/' + prefix + '.dyn_q0'
            ))
            tmp_path = os.path.join(self._FOLDER_DYNAMICAL_MATRIX, 'dynamical-matrix-' + label)
            remote_copy_list.append((
                parent_folder_ph.computer.uuid, os.path.join(parent_folder_ph.get_remote_path(),
                                                             tmp_path), 'save/' + prefix + '.dyn_q' + label
            ))

            if iqpt == 1:
                tmp_path = os.path.join(self._OUTPUT_SUBFOLDER, '_ph0/' + prefix + '.dvscf*')
                remote_copy_list.append((
                    parent_folder_ph.computer.uuid, os.path.join(parent_folder_ph.get_remote_path(),
                                                                 tmp_path), 'save/' + prefix + '.dvscf_q' + label
                ))
                tmp_path = os.path.join(self._OUTPUT_SUBFOLDER, '_ph0/' + prefix + '.phsave')
                remote_copy_list.append((
                    parent_folder_ph.computer.uuid, os.path.join(parent_folder_ph.get_remote_path(), tmp_path), 'save/'
                ))
            else:
                tmp_path = os.path.join(
                    self._OUTPUT_SUBFOLDER, '_ph0/' + prefix + '.q_' + label + '/' + prefix + '.dvscf*'
                )
                remote_copy_list.append((
                    parent_folder_ph.computer.uuid, os.path.join(parent_folder_ph.get_remote_path(),
                                                                 tmp_path), 'save/' + prefix + '.dvscf_q' + label
                ))

        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = (list(settings.pop('CMDLINE', [])) + ['-in', self.metadata.options.input_filename])
        codeinfo.stdout_name = self.metadata.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid

        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.remote_symlink_list = remote_symlink_list

        # Retrieve by default the output file
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(self.metadata.options.output_filename)
        calcinfo.retrieve_list += settings.pop('ADDITIONAL_RETRIEVE_LIST', [])

        if settings:
            unknown_keys = ', '.join(list(settings.keys()))
            raise exceptions.InputValidationError(f'`settings` contained unexpected keys: {unknown_keys}')

        return calcinfo
