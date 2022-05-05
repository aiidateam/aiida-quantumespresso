# -*- coding: utf-8 -*-
"""Calculation function to compute a k-point mesh for a structure with a guaranteed minimum k-point distance."""
# pylint: disable=unused-import
import warnings

from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import create_kpoints_from_distance

warnings.warn(
    'This module is deprecated and will be removed soon.\nPlease use instead the new module:\n'
    'from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import create_kpoints_from_distance'
    "\nOr use the entry point with the factory: CalculationFactory('quantumespresso.create_kpoints_from_distance')",
    FutureWarning
)
