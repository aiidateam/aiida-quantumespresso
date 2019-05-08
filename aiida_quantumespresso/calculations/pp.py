# -*- coding: utf-8 -*-

from __future__ import absolute_import
from aiida.orm import RemoteData, FolderData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class PpCalculation(NamelistsCalculation):
    """
    `pp.x` code of the Quantum ESPRESSO distribution, handles the
    post-processing of charge-densities, potentials, ...
    For more information, refer to http://www.quantum-espresso.org/
    """

    _FILPLOT = "aiida.filplot"

    _default_namelists = ["INPUTPP", "PLOT"]

    _blocked_keywords = [
        ("INPUTPP", "outdir", NamelistsCalculation._OUTPUT_SUBFOLDER),
        ("INPUTPP", "prefix", NamelistsCalculation._PREFIX),
        ("INPUTPP", "filplot", NamelistsCalculation._FILPLOT),
    ]

    _default_parser = None

    _internal_retrieve_list = [_FILPLOT]

    @classmethod
    def define(cls, spec):
        super(PpCalculation, cls).define(spec)
        spec.input(
            "parent_folder",
            valid_type=(RemoteData, FolderData),
            required=True,
            help="Output folder of a PW calculation",
        )
