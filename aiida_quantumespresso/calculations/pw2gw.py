# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the pw2gw.x code of Quantum ESPRESSO."""
from __future__ import absolute_import

import os

from aiida import orm
from aiida.common import datastructures
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class Pw2gwCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the pw2gw.x code of Quantum ESPRESSO."""

    _INPUT_PSEUDOFOLDER = './pseudo/'

    _EPS_X = 'epsX.dat'
    _EPS_Y = 'epsY.dat'
    _EPS_Z = 'epsZ.dat'
    _EPS_TOT = 'epsTOT.dat'

    _default_namelists = ['INPUTPP']
    _internal_retrieve_list = [_EPS_X, _EPS_Y, _EPS_Z, _EPS_TOT]
    _blocked_keywords = [
        ('INPUTPP', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),
        ('INPUTPP', 'prefix', NamelistsCalculation._PREFIX),
    ]

    @classmethod
    def define(cls, spec):
        # yapf: disable
        super(Pw2gwCalculation, cls).define(spec)
        spec.input('parent_folder', valid_type=orm.RemoteData,
            help='Output folder of a completed `PwCalculation`')

        spec.output('output_parameters', valid_type=orm.Dict,
            help='The `output_parameters` output node of the successful calculation.`')
        spec.output('eps', valid_type=orm.ArrayData,
            help='The `eps` output node containing 5 arrays `energy`, `epsX`, `epsY`, `epsZ`, `epsTOT`'
            )

        spec.exit_code(300, 'ERROR_NO_RETRIEVED_FOLDER',
            message='The retrieved folder data node could not be accessed.')
        spec.exit_code(302, 'ERROR_OUTPUT_STDOUT_MISSING',
            message='The retrieved folder did not contain the required stdout output file.')
        spec.exit_code(305, 'ERROR_OUTPUT_FILES',
            message='The eps*.dat output files could not be read or parsed.')
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(311, 'ERROR_OUTPUT_STDOUT_PARSE',
            message='The stdout output file could not be parsed.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete probably because the calculation got interrupted.')
        spec.exit_code(330, 'ERROR_OUTPUT_FILES_INVALID_FORMAT',
            message='The eps*.dat output files do not have the expected shape (N, 2).')
        spec.exit_code(331, 'ERROR_OUTPUT_FILES_ENERGY_MISMATCH',
            message='The eps*.dat output files contains different values of energies.')
        spec.exit_code(350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION',
            message='The parser raised an unexpected exception.')

    def prepare_for_submission(self, folder):
        # yapf: disable
        calcinfo = super(Pw2gwCalculation, self).prepare_for_submission(folder)

        calcinfo.codes_run_mode = datastructures.CodeRunMode.SERIAL

        parent_calc_folder = self.inputs.parent_folder
        calcinfo.remote_copy_list.append((
            parent_calc_folder.computer.uuid,
            os.path.join(parent_calc_folder.get_remote_path(), self._INPUT_PSEUDOFOLDER),
            self._INPUT_PSEUDOFOLDER
        ))

        return calcinfo
