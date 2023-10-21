# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the q2r.x code of Quantum ESPRESSO."""

from aiida import orm

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
