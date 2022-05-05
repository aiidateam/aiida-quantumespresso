# -*- coding: utf-8 -*-
"""Calcfunction to primitivize a structure and return high symmetry k-point path through its Brillouin zone."""
# pylint: disable=unused-import
import warnings

from aiida_quantumespresso.calculations.functions.seekpath_structure_analysis import seekpath_structure_analysis

warnings.warn(
    'This module is deprecated and will be removed soon.\nPlease use instead the new module:\n'
    'from aiida_quantumespresso.calculations.functions.seekpath_structure_analysis import seekpath_structure_analysis'
    "\nOr use the entry point with the factory: CalculationFactory('quantumespresso.seekpath_structure_analysis')",
    FutureWarning
)
