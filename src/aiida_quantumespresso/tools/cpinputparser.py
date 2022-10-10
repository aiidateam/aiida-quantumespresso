# -*- coding: utf-8 -*-
"""Utilities to parse Quantum ESPRESSO cp.x input files into AiiDA nodes or builders."""
from qe_tools.parsers import CpInputFile as BaseCpInputFile

from .base import StructureParseMixin


class CpInputFile(StructureParseMixin, BaseCpInputFile):  # pylint: disable=too-few-public-methods
    """Parser of Quantum ESPRESSO cp.x input file into AiiDA nodes.

    .. note:: This mixes in :class:`~aiida_quantumespresso.tools.base.StructureParseMixin` which adds the functionality
        to parse a :class:`~aiida.orm.nodes.data.structure.StructureData` from the input file, instead of a plain
        dictionary returned by ``qe_tools.parsers.qeinputparser.get_structure_from_qeinput``. Note that one cannot
        directly add this functionality to a sub class of ``~qe_tools.parsers.qeinputparser.QeInputFile`` and then
        subsequently sub class that here, because the ``~qe_tools.parsers.qeinputparser.CpInputFile`` is also
        required and sub classing both leads to problems with the MRO.
    """
