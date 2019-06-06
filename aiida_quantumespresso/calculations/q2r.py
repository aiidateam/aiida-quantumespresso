# -*- coding: utf-8 -*-
from __future__ import absolute_import

import os

from aiida import orm
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation
from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.data.forceconstants import ForceconstantsData


class Q2rCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the q2r.x code of Quantum ESPRESSO."""

    _FORCE_CONSTANTS_NAME = 'real_space_force_constants.dat'
    _OUTPUT_SUBFOLDER = PhCalculation._FOLDER_DYNAMICAL_MATRIX
    _INPUT_SUBFOLDER = os.path.join('.', PhCalculation._FOLDER_DYNAMICAL_MATRIX)

    _default_namelists = ['INPUT']
    _blocked_keywords = [
        ('INPUT', 'fildyn', PhCalculation._OUTPUT_DYNAMICAL_MATRIX_PREFIX),
        ('INPUT', 'flfrc', _FORCE_CONSTANTS_NAME),
    ]

    _internal_retrieve_list = [_FORCE_CONSTANTS_NAME]
    _default_parser = 'quantumespresso.q2r'

    @classmethod
    def define(cls, spec):
        super(Q2rCalculation, cls).define(spec)
        spec.input('parent_folder', valid_type=(orm.RemoteData, orm.FolderData), required=True)
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('forceconstants', valid_type=ForceconstantsData)
        spec.default_output_node = 'output_parameters'
        spec.exit_code(
            100, 'ERROR_NO_RETRIEVED_FOLDER', message='The retrieved folder data node could not be accessed.')
        spec.exit_code(
            110, 'ERROR_READING_OUTPUT_FILE', message='The output file could not be read from the retrieved folder.')
        spec.exit_code(
            130, 'ERROR_JOB_NOT_DONE', message='The computation did not finish properly ("JOB DONE" not found).')
