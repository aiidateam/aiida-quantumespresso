# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the pw2gw.x code of Quantum ESPRESSO."""
from __future__ import absolute_import

from aiida.orm import RemoteData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class Pw2gwCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the pw2gw.x code of Quantum ESPRESSO."""

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
        spec.input('parent_folder', valid_type=RemoteData,
            help='Output folder of a completed `PwCalculation`')

        spec.exit_code(200, 'ERROR_NO_RETRIEVED_FOLDER',
            message='The retrieved folder data node could not be accessed.')
        spec.exit_code(300, 'ERROR_OUTPUT_FILES',
            message='The eps*.dat output files could not be read or parsed.')
