# -*- coding: utf-8 -*-
"""Tests for the utility functions for the protocols."""

import pytest

from aiida_quantumespresso.common.types import SpinType


def test_recursive_merge():
    """Test the recursive merge function."""
    from aiida_quantumespresso.workflows.protocols.utils import recursive_merge

    left_dict = {'a': {'b': 1, 'c': {'d': 2}}, 'g': 3}
    right_dict = {'a': {'c': {'d': 'D'}}, 'e': {'f': 'F'}}
    merged = recursive_merge(left_dict, right_dict)

    assert right_dict == {'a': {'c': {'d': 'D'}}, 'e': {'f': 'F'}}

    assert merged == {'a': {'b': 1, 'c': {'d': 'D'}}, 'e': {'f': 'F'}, 'g': 3}


# Note: for the `pseudo_family` text fixture, the z_valence is 4 for every element.
# This can lead to some nonsensical magnetizations, but that is fine for testing
# purposes.
@pytest.mark.parametrize(
    'structure_id,initial_magnetic_moments,spin_type,expected_magnetization',
    (
        (
            'silicon',
            None,
            SpinType.COLLINEAR,
            {
                'starting_magnetization': {
                    'Si': 0.1
                },
                'angle1': None,
                'angle2': None
            },
        ),
        (
            'cobalt-prim',
            None,
            SpinType.COLLINEAR,
            {
                'starting_magnetization': {
                    'Co': 1.25
                },
                'angle1': None,
                'angle2': None
            },
        ),
        (
            'cobalt-prim',
            None,
            SpinType.NON_COLLINEAR,
            {
                'starting_magnetization': {
                    'Co': 1.25
                },
                'angle1': {
                    'Co': 0
                },
                'angle2': {
                    'Co': 0
                }
            },
        ),
        (
            'cobalt-prim',
            None,
            SpinType.SPIN_ORBIT,
            {
                'starting_magnetization': {
                    'Co': 1.25
                },
                'angle1': {
                    'Co': 0
                },
                'angle2': {
                    'Co': 0
                }
            },
        ),
        (
            'cobalt-prim',
            {
                'Co': 3
            },
            SpinType.COLLINEAR,
            {
                'starting_magnetization': {
                    'Co': 0.75
                },
                'angle1': None,
                'angle2': None
            },
        ),
        (
            'cobalt-prim',
            {
                'Co': (1, 2, 3)
            },
            SpinType.NON_COLLINEAR,
            {
                'starting_magnetization': {
                    'Co': 0.25
                },
                'angle1': {
                    'Co': 2
                },
                'angle2': {
                    'Co': 3
                }
            },
        ),
    ),
)
def test_get_magnetization(
    generate_structure,
    structure_id,
    initial_magnetic_moments,
    spin_type,
    expected_magnetization,
):
    """Test the `get_magnetization` function."""
    from aiida_quantumespresso.workflows.protocols.utils import get_magnetization
    structure = generate_structure(structure_id)
    z_valences = {kind: 4.0 for kind in structure.get_kind_names()}

    magnetization = get_magnetization(structure, z_valences, initial_magnetic_moments, spin_type)

    assert magnetization == expected_magnetization


@pytest.mark.parametrize(
    'structure_id,z_valences,initial_magnetic_moments,spin_type,expected_error,error_message',
    (
        ('silicon', {}, {
            'Si': 1.0
        }, SpinType.COLLINEAR, ValueError, '`z_valences` needs one value for each of'),
        (
            'silicon', {
                'Si': 4.0
            }, {}, SpinType.COLLINEAR, ValueError, '`initial_magnetic_moments` needs one value for each of'
        ),
        ('silicon', {
            'Si': 4.0
        }, {
            'Si': (1, 2, 3)
        }, SpinType.COLLINEAR, TypeError, 'Spin type is set to '),
        (
            'silicon', {
                'Si': 4.0
            }, {
                'Si': 'zero'
            }, SpinType.COLLINEAR, TypeError, 'Unrecognised type for magnetic moment'
        ),
    ),
)
def test_get_magnetization_failure(
    generate_structure, structure_id, z_valences, initial_magnetic_moments, spin_type, expected_error, error_message
):
    """Test the `get_magnetization` function."""
    from aiida_quantumespresso.workflows.protocols.utils import get_magnetization
    with pytest.raises(expected_error, match=error_message):
        get_magnetization(generate_structure(structure_id), z_valences, initial_magnetic_moments, spin_type)
