# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the matdyn.x code of Quantum ESPRESSO."""
from pathlib import Path
import warnings

from aiida import orm

from aiida_quantumespresso.calculations import _uppercase_dict
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation
from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.data.force_constants import ForceConstantsData


class MatdynCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the matdyn.x code of Quantum ESPRESSO."""

    _PHONON_FREQUENCIES_NAME = 'phonon_frequencies.dat'
    _PHONON_MODES_NAME = 'phonon_displacements.dat'
    _PHONON_DOS_NAME = 'phonon_dos.dat'

    _default_namelists = ['INPUT']
    _blocked_keywords = [
        ('INPUT', 'flfrq', _PHONON_FREQUENCIES_NAME),  # output frequencies
        ('INPUT', 'flvec', _PHONON_MODES_NAME),  # output displacements
        ('INPUT', 'fldos', _PHONON_DOS_NAME),  # output density of states
    ]

    _internal_retrieve_list = [_PHONON_FREQUENCIES_NAME, _PHONON_DOS_NAME]
    _default_parser = 'quantumespresso.matdyn'

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('force_constants', valid_type=ForceConstantsData, required=True)
        spec.input('kpoints', valid_type=orm.KpointsData, help='Kpoints on which to calculate the phonon frequencies.')
        spec.input('parent_folder', valid_type=orm.RemoteData, required=False)
        spec.inputs.validator = cls._validate_inputs

        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('output_phonon_bands', valid_type=orm.BandsData, required=False)
        spec.output('output_phonon_dos', valid_type=orm.XyData, required=False)
        spec.default_output_node = 'output_parameters'

        spec.exit_code(330, 'ERROR_OUTPUT_FREQUENCIES',
            message='The output frequencies file could not be read from the retrieved folder.')
        spec.exit_code(334, 'ERROR_OUTPUT_DOS',
            message='The output DOS file could not be read from the retrieved folder.')
        spec.exit_code(410, 'ERROR_OUTPUT_KPOINTS_MISSING',
            message='Number of kpoints not found in the output data')
        spec.exit_code(411, 'ERROR_OUTPUT_KPOINTS_INCOMMENSURATE',
            message='Number of kpoints in the inputs is not commensurate with those in the output')
        # yapf: enable

    @staticmethod
    def _validate_inputs(value, _):
        """Validate the top level namespace."""
        if 'parameters' in value:
            parameters = value['parameters'].get_dict()
        else:
            parameters = {'INPUT': {}}

        if 'INPUT' not in parameters:
            return 'Required namelist `INPUT` not in `parameters` input.'

        if parameters['INPUT'].get('flfrc', None) is not None:
            warnings.warn(
                '`INPUT.flfrc` is set automatically when the `force_constants` input is specified.'
                'There is no need to specify this input, and its value will be overridden.'
            )
        if parameters['INPUT'].get('q_in_cryst_coord', None) is not None:
            warnings.warn(
                '`INPUT.q_in_cryst_coords` is always set to `.true.` if `kpoints` input corresponds to list.'
                'There is no need to specify this input, and its value will be overridden.'
            )

        if 'parent_folder' in value and not parameters['INPUT'].get('la2F', False):
            return (
                'The `parent_folder` input is only used to calculate the el-ph coefficients but `la2F` is not set '
                'to `.true.` in input `parameters`'
            )

    def generate_input_file(self, parameters):  # pylint: disable=arguments-differ
        """Generate namelist input_file content given a dict of parameters.

        :param parameters: 'dict' containing the fortran namelists and parameters to be used.
          e.g.: {'CONTROL':{'calculation':'scf'}, 'SYSTEM':{'ecutwfc':30}}
        :return: 'str' containing the input_file content a plain text.
        """
        kpoints = self.inputs.kpoints
        append_string = ''

        parameters.setdefault('INPUT', {})['flfrc'] = self.inputs.force_constants.filename

        # Calculating DOS requires (nk1,nk2,nk3), see
        # https://gitlab.com/QEF/q-e/-/blob/b231a0d0174ad1853f191160389029aa14fba6e9/PHonon/PH/matdyn.f90#L82
        if parameters['INPUT'].get('dos', False):
            kpoints_mesh = kpoints.get_kpoints_mesh()[0]
            parameters['INPUT']['nk1'] = kpoints_mesh[0]
            parameters['INPUT']['nk2'] = kpoints_mesh[1]
            parameters['INPUT']['nk3'] = kpoints_mesh[2]
        else:
            parameters['INPUT']['q_in_cryst_coord'] = True
            try:
                kpoints_list = kpoints.get_kpoints()
            except AttributeError:
                kpoints_list = kpoints.get_kpoints_mesh(print_list=True)

            kpoints_string = [f'{len(kpoints_list)}']
            for kpoint in kpoints_list:
                kpoints_string.append('{:18.10f} {:18.10f} {:18.10f}'.format(*kpoint))  # pylint: disable=consider-using-f-string
            append_string = '\n'.join(kpoints_string) + '\n'

        file_content = super().generate_input_file(parameters)

        return file_content + append_string

    def prepare_for_submission(self, folder):
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        After calling the method of the parent `NamelistsCalculation` class, the input parameters are checked to see
        if the `la2F` tag is set to true. In this case the remote symlink or copy list is set to the electron-phonon
        directory, depending on the settings.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :class:`~aiida.common.datastructures.CalcInfo` instance.
        """
        calcinfo = super().prepare_for_submission(folder)

        force_constants = self.inputs.force_constants
        calcinfo.local_copy_list.append((force_constants.uuid, force_constants.filename, force_constants.filename))

        if 'settings' in self.inputs:
            settings = _uppercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
        else:
            settings = {}

        if 'parameters' in self.inputs:
            parameters = _uppercase_dict(self.inputs.parameters.get_dict(), dict_name='parameters')
        else:
            parameters = {'INPUT': {}}

        source = self.inputs.get('parent_folder', None)

        if source is not None and parameters['INPUT'].get('la2F', False):

            # pylint: disable=protected-access
            dirpath = Path(source.get_remote_path()) / PhCalculation._FOLDER_ELECTRON_PHONON
            remote_list = [(source.computer.uuid, str(dirpath), PhCalculation._FOLDER_ELECTRON_PHONON)]

            # For el-ph calculations, _only_ the `elph_dir` should be copied from the parent folder
            if settings.pop('PARENT_FOLDER_SYMLINK', False):
                calcinfo.remote_symlink_list = remote_list
            else:
                calcinfo.remote_copy_list = remote_list

            calcinfo.retrieve_list += [f'a2F.dos{i}' for i in range(1, 11)]
            calcinfo.retrieve_list.append('lambda')

        return calcinfo
