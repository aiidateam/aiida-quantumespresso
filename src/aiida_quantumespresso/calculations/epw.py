# -*- coding: utf-8 -*-
"""Plugin to create a Quantum Espresso epw.x input file."""
from pathlib import Path

from aiida import orm
from aiida.common import datastructures, exceptions

from aiida_quantumespresso.calculations import _lowercase_dict, _uppercase_dict
from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.calculations.pw import PwCalculation
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
    _FOLDER_SAVE = 'save'
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
        spec.input('parent_folder_nscf', required=False, valid_type=orm.RemoteData,
                   help='the folder of a completed nscf `PwCalculation`')
        spec.input('parent_folder_ph', required=False, valid_type=orm.RemoteData,
                   help='the folder of a completed `PhCalculation`')
        spec.input('parent_folder_epw', required=False, valid_type=(orm.RemoteData, orm.RemoteStashFolderData),
                   help='folder that contains all files required to restart an `EpwCalculation`')
        # yapf: enable

    def prepare_for_submission(self, folder):
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :class:`~aiida.common.datastructures.CalcInfo` instance.
        """

        # pylint: disable=too-many-statements,too-many-branches, protected-access

        def test_offset(offset):
            """Check if the grid has an offset."""
            if any(i != 0. for i in offset):
                raise NotImplementedError(
                    'Computation of electron-phonon on a mesh with non zero offset is not implemented, '
                    'at the level of epw.x'
                )

        local_copy_list = []
        remote_copy_list = []
        remote_symlink_list = []

        parameters = _uppercase_dict(self.inputs.parameters.get_dict(), dict_name='parameters')
        parameters = {k: _lowercase_dict(v, dict_name=k) for k, v in parameters.items()}

        if 'INPUTEPW' not in parameters:
            raise exceptions.InputValidationError('required namelist INPUTEPW not specified')

        if 'settings' in self.inputs:
            settings = _uppercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
        else:
            settings = {}

        remote_list = remote_symlink_list if settings.pop(
            'PARENT_FOLDER_SYMLINK', self._default_symlink_usage
        ) else remote_copy_list

        if 'parent_folder_nscf' in self.inputs:
            parent_folder_nscf = self.inputs.parent_folder_nscf

            remote_list.append((
                parent_folder_nscf.computer.uuid,
                Path(parent_folder_nscf.get_remote_path(), PwCalculation._OUTPUT_SUBFOLDER).as_posix(),
                self._OUTPUT_SUBFOLDER,
            ))

        if 'parent_folder_ph' in self.inputs:
            parent_folder_ph = self.inputs.parent_folder_ph

            # Create the save folder with dvscf and dyn files
            folder.get_subfolder(self._FOLDER_SAVE, create=True)

            # List of IBZ q-point to be added below EPW. To be removed when removed from EPW.
            qibz_ar = []
            for key, value in sorted(parent_folder_ph.creator.outputs.output_parameters.get_dict().items()):
                if key.startswith('dynamical_matrix_'):
                    qibz_ar.append(value['q_point'])

            nqpt = len(qibz_ar)

            # Append the required contents of the `save` folder to the remove copy list, copied from the `ph.x`
            # calculation

            prefix = self._PREFIX
            outdir = self._OUTPUT_SUBFOLDER
            fildvscf = PhCalculation._DVSCF_PREFIX
            fildyn = PhCalculation._OUTPUT_DYNAMICAL_MATRIX_PREFIX

            ph_path = Path(parent_folder_ph.get_remote_path())

            remote_list.append(
                (parent_folder_ph.computer.uuid, Path(ph_path, outdir, '_ph0', f'{prefix}.phsave').as_posix(), 'save')
            )

            for iqpt in range(1, nqpt + 1):
                remote_list.append((
                    parent_folder_ph.computer.uuid,
                    Path(ph_path, outdir, '_ph0', '' if iqpt == 1 else f'{prefix}.q_{iqpt}',
                         f'{prefix}.{fildvscf}1').as_posix(), Path('save', f'{prefix}.dvscf_q{iqpt}').as_posix()
                ))
                # remote_copy_list.append((
                #     parent_folder_ph.computer.uuid,
                #     Path(
                #     ph_path, outdir, '_ph0', '' if iqpt == 1 else f'{prefix}.q_{iqpt}', f'{prefix}.{fildvscf}_paw1'
                #     ).as_posix(),
                #     Path('save', f"{prefix}.dvscf_paw_q{iqpt}").as_posix()
                # ))
                remote_list.append((
                    parent_folder_ph.computer.uuid, Path(ph_path, f'{fildyn}{iqpt}').as_posix(),
                    Path('save', f'{prefix}.dyn_q{iqpt}').as_posix()
                ))

        if 'parent_folder_epw' in self.inputs:

            parent_folder_epw = self.inputs.parent_folder_epw
            if isinstance(parent_folder_epw, orm.RemoteStashFolderData):
                epw_path = Path(parent_folder_epw.target_basepath)
            else:
                epw_path = Path(parent_folder_epw.get_remote_path())

            vme_fmt_dict = {
                'dipole': 'dmedata.fmt',
                'wannier': 'vmedata.fmt',
            }

            for filename in (
                'crystal.fmt', 'epwdata.fmt', vme_fmt_dict[parameters['INPUTEPW']['vme']], f'{self._PREFIX}.kgmap',
                f'{self._PREFIX}.kmap', f'{self._PREFIX}.ukk', self._OUTPUT_SUBFOLDER, self._FOLDER_SAVE
            ):
                remote_list.append(
                    (parent_folder_epw.computer.uuid, Path(epw_path, filename).as_posix(), Path(filename).as_posix())
                )

        parameters['INPUTEPW']['outdir'] = self._OUTPUT_SUBFOLDER
        parameters['INPUTEPW']['dvscf_dir'] = self._FOLDER_SAVE
        parameters['INPUTEPW']['prefix'] = self._PREFIX

        try:
            mesh, offset = self.inputs.qpoints.get_kpoints_mesh()
            test_offset(offset)
            parameters['INPUTEPW']['nq1'] = mesh[0]
            parameters['INPUTEPW']['nq2'] = mesh[1]
            parameters['INPUTEPW']['nq3'] = mesh[2]
        except NotImplementedError as exception:
            raise exceptions.InputValidationError('Cannot get the coarse q-point grid') from exception

        try:
            mesh, offset = self.inputs.kpoints.get_kpoints_mesh()
            test_offset(offset)
            parameters['INPUTEPW']['nk1'] = mesh[0]
            parameters['INPUTEPW']['nk2'] = mesh[1]
            parameters['INPUTEPW']['nk3'] = mesh[2]
        except NotImplementedError as exception:
            raise exceptions.InputValidationError('Cannot get the coarse k-point grid') from exception

        try:
            mesh, offset = self.inputs.qfpoints.get_kpoints_mesh()
            test_offset(offset)
            parameters['INPUTEPW']['nqf1'] = mesh[0]
            parameters['INPUTEPW']['nqf2'] = mesh[1]
            parameters['INPUTEPW']['nqf3'] = mesh[2]
        except NotImplementedError as exception:
            raise exceptions.InputValidationError('Cannot get the fine q-point grid') from exception

        try:
            mesh, offset = self.inputs.kfpoints.get_kpoints_mesh()
            test_offset(offset)
            parameters['INPUTEPW']['nkf1'] = mesh[0]
            parameters['INPUTEPW']['nkf2'] = mesh[1]
            parameters['INPUTEPW']['nkf3'] = mesh[2]
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

        with folder.open(self.metadata.options.input_filename, 'w') as infile:
            for namelist_name in namelists_toprint:
                infile.write(f'&{namelist_name}\n')
                # namelist content; set to {} if not present, so that we leave an empty namelist
                namelist = parameters.pop(namelist_name, {})
                for key, value in sorted(namelist.items()):
                    infile.write(convert_input_to_namelist_entry(key, value))
                infile.write('/\n')

        if parameters:
            raise exceptions.InputValidationError(
                'The following namelists are specified in parameters, but are not valid namelists for the current type '
                f'of calculation: {",".join(list(parameters.keys()))}'
            )

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
