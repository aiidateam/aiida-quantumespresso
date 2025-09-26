# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the q2r.x code of Quantum ESPRESSO."""
from pathlib import Path

from aiida import orm

from aiida_quantumespresso.calculations import _uppercase_dict
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation
from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.data.force_constants import ForceConstantsData


class Q2rCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the q2r.x code of Quantum ESPRESSO."""

    _FORCE_CONSTANTS_NAME = 'real_space_force_constants.dat'
    _OUTPUT_SUBFOLDER = PhCalculation._FOLDER_DYNAMICAL_MATRIX  # pylint: disable=protected-access
    _INPUT_SUBFOLDER = PhCalculation._FOLDER_DYNAMICAL_MATRIX  # pylint: disable=protected-access
    _default_parent_output_folder = PhCalculation._FOLDER_DYNAMICAL_MATRIX  # pylint: disable=protected-access

    _default_namelists = ['INPUT']
    _blocked_keywords = [
        ('INPUT', 'fildyn', PhCalculation._OUTPUT_DYNAMICAL_MATRIX_PREFIX),  # pylint: disable=protected-access
        ('INPUT', 'flfrc', _FORCE_CONSTANTS_NAME),  # Real space force constants
    ]

    _internal_retrieve_list = [_FORCE_CONSTANTS_NAME]
    _default_parser = 'quantumespresso.q2r'

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('parent_folder', valid_type=(orm.RemoteData, orm.FolderData), required=True)
        spec.output('force_constants', valid_type=ForceConstantsData)
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.exit_code(330, 'ERROR_READING_FORCE_CONSTANTS_FILE',
            message='The force constants file could not be read.')
        # yapf: enable

    def prepare_for_submission(self, folder):
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        After calling the method of the parent `NamelistsCalculation` class, the input parameters are checked to see
        if the `la2F` tag is set to true. In this case the electron-phonon directory is added to the remote symlink or
        copy list, depending on the settings.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :class:`~aiida.common.datastructures.CalcInfo` instance.
        """
        calcinfo = super().prepare_for_submission(folder)

        if 'settings' in self.inputs:
            settings = _uppercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
        else:
            settings = {}

        if 'parameters' in self.inputs:
            parameters = _uppercase_dict(self.inputs.parameters.get_dict(), dict_name='parameters')
        else:
            parameters = {'INPUT': {}}

        parent_folder = self.inputs.get('parent_folder', None)

        if parent_folder is not None and 'parameters' in self.inputs:

            if parameters.get('INPUT').get('la2F', False):

                symlink = settings.pop('PARENT_FOLDER_SYMLINK', False)
                remote_list = calcinfo.remote_symlink_list if symlink else calcinfo.remote_copy_list

                # pylint: disable=protected-access
                dirpath = Path(parent_folder.get_remote_path()) / PhCalculation._FOLDER_ELECTRON_PHONON
                remote_list.append((parent_folder.computer.uuid, str(dirpath), PhCalculation._FOLDER_ELECTRON_PHONON))

        return calcinfo
