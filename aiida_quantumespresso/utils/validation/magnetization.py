# -*- coding: utf-8 -*-
"""Utilities for validating various data structures."""


def validate_parent_calculation(calculation):
    """Validate whether the given calculation is a valid parent calculation for a HpCalculation.

    :param calculation: the calculation to validate
    :raises ValueError: if the calculation is not a valid parent calculation for a HpCalculation
    """
    from aiida.plugins import CalculationFactory

    PwCalculation = CalculationFactory('quantumespresso.pw')

    if not hasattr(calculation, 'process_class') or calculation.process_class is not PwCalculation:
        raise ValueError(f'parent calculation is not of type `PwCalculation` but {calculation}')

    try:
        parameters = calculation.inputs.parameters.get_dict()
    except AttributeError as exception:
        raise ValueError('could not retrieve the input parameters node') from exception

    lda_plus_u = parameters.get('SYSTEM', {}).get('lda_plus_u', None)
    hubbard_u = parameters.get('SYSTEM', {}).get('hubbard_u', {})
    hubbard_v = parameters.get('SYSTEM', {}).get('hubbard_v', {})
    hubbard_parameters = parameters.get('SYSTEM', {}).get('hubbard_parameters', None)

    if lda_plus_u is not True:
        raise ValueError('the parent calculation did not set `lda_plus_u=True`')

    if not hubbard_u and not hubbard_v and not hubbard_parameters:
        raise ValueError('the parent calculation did not specify any Hubbard U or V parameters')

    try:
        structure = calculation.inputs.structure
    except AttributeError as exception:
        raise ValueError('could not retrieve the input structure node') from exception

    validate_structure_kind_order(structure, list(hubbard_u.keys()))


def validate_structure_kind_order(structure, hubbard_kinds):
    """Determine whether the kinds in the structure node have the right order for the given list of Hubbard U kinds.

    For the order to be right, means for the Hubbard kinds to come first in the list of kinds of the structure.

    :param structure: StructureData node
    :param hubbard_kinds: a list of Hubbard kinds
    """
    for kind in structure.kinds:

        if not hubbard_kinds:
            return

        if kind.name in hubbard_kinds:
            hubbard_kinds.remove(kind.name)
        else:
            raise ValueError('the structure does not have the right kind order')
