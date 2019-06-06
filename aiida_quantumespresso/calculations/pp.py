# -*- coding: utf-8 -*-
from __future__ import absolute_import

from aiida.orm import RemoteData, FolderData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class PpCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the pp.x code of Quantum ESPRESSO."""

    _FILPLOT = 'aiida.filplot'

    _default_namelists = ['INPUTPP', 'PLOT']
    _internal_retrieve_list = [_FILPLOT]
    _blocked_keywords = [
        ('INPUTPP', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),
        ('INPUTPP', 'prefix', NamelistsCalculation._PREFIX),
        ('INPUTPP', 'filplot', _FILPLOT),
    ]

    @classmethod
    def define(cls, spec):
        super(PpCalculation, cls).define(spec)
        spec.input('parent_folder', valid_type=(RemoteData, FolderData),
            help='Output folder of a completed `PwCalculation`')
