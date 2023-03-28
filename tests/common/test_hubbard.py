# -*- coding: utf-8 -*-
"""Tests for :py:mod:`~aiida_quantumespresso.common.hubbard`."""
# pylint: disable=redefined-outer-name
from copy import deepcopy

from pydantic import ValidationError
import pytest

from aiida_quantumespresso.common.hubbard import Hubbard, HubbardParameters

VALID_PARAMETERS = {
    'atom_index': 0,
    'atom_manifold': '3d',
    'neighbour_index': 1,
    'neighbour_manifold': '2p',
    'translation': (0, 0, 0),
    'value': 5.0,
    'hubbard_type': 'U',
}


@pytest.fixture
def get_hubbard_parameters():
    """Return an `HubbardParameters` intstance."""

    def _get_hubbard_parameters(overrides=None):
        """Return an `HubbardParameters` intstance."""
        inputs = deepcopy(VALID_PARAMETERS)

        if overrides:
            inputs.update(overrides)

        return HubbardParameters(**inputs)

    return _get_hubbard_parameters


@pytest.fixture
def get_hubbard():
    """Return an `Hubbard` intstance."""

    def _get_hubbard():
        """Return an `Hubbard` intstance."""
        param = HubbardParameters(**VALID_PARAMETERS)

        return Hubbard(parameters=[param, param])

    return _get_hubbard


def test_safe_hubbard_parameters(get_hubbard_parameters):
    """Test valid inputs are stored correctly for py:meth:`HubbardParameters`."""
    params = get_hubbard_parameters().dict()
    assert params == VALID_PARAMETERS


def test_from_to_list_parameters(get_hubbard_parameters):
    """Test py:meth:`HubbardParameters.to_tuple` and py:meth:`HubbardParameters.from_tuple`."""
    param = get_hubbard_parameters()
    hp_tuple = (0, '3d', 1, '2p', 5.0, (0, 0, 0), 'U')
    assert param.to_tuple() == hp_tuple
    param = HubbardParameters.from_tuple(hp_tuple)
    assert param.dict() == VALID_PARAMETERS


@pytest.mark.parametrize(
    'overrides', [{
        'atom_index': 0
    }, {
        'atom_manifold': '3d-2p'
    }, {
        'translation': (0, -1, +1)
    }, {
        'hubbard_type': 'B'
    }]
)
def test_valid_hubbard_parameters(get_hubbard_parameters, overrides):
    """Test valid inputs for py:meth:`HubbardParameters`."""
    hp_dict = get_hubbard_parameters(overrides=overrides).dict()
    new_dict = deepcopy(VALID_PARAMETERS)
    new_dict.update(overrides)
    assert hp_dict == new_dict


@pytest.mark.parametrize(('overrides', 'match'), (
    ({
        'atom_index': -1
    }, r'ensure this value is greater than or equal to 0'),
    (
        {
            'atom_index': 0.5
        },
        r'value is not a valid integer',
    ),
    (
        {
            'atom_manifold': '3z'
        },
        r'invalid manifold symbol z',
    ),
    (
        {
            'atom_manifold': '3d2p'
        },
        r'invalid length ``4``. Only 2 or 5',
    ),
    (
        {
            'atom_manifold': '3d-3p-2s'
        },
        r'ensure this value has at most 5 characters',
    ),
    (
        {
            'translation': (0, 0)
        },
        r'wrong tuple length 2, expected 3',
    ),
    (
        {
            'translation': (0, 0, 0, 0)
        },
        r'wrong tuple length 4, expected 3',
    ),
    (
        {
            'translation': (0, 0, -1.5)
        },
        r'value is not a valid integer',
    ),
    (
        {
            'hubbard_type': 'L'
        },
        r"permitted: 'Ueff', 'U', 'V', 'J', 'B', 'E2', 'E3'",
    ),
))
def test_invalid_hubbard_parameters(get_hubbard_parameters, overrides, match):
    """Test invalid inputs for py:meth:`HubbardParameters`."""
    with pytest.raises(ValidationError, match=match):
        get_hubbard_parameters(overrides=overrides)


def test_from_to_list_hubbard(get_hubbard):
    """Test py:meth:`Hubbard.to_list` and py:meth:`Hubbard.from_list`."""
    hubbard = get_hubbard()
    hp_tuple = (0, '3d', 1, '2p', 5.0, (0, 0, 0), 'U')

    hubbard_list = [hp_tuple, hp_tuple]
    assert hubbard.to_list() == hubbard_list

    hubbard = Hubbard.from_list(hubbard_list)
    assert hubbard.dict() == {
        'parameters': [VALID_PARAMETERS, VALID_PARAMETERS],
        'projectors': 'ortho-atomic',
        'formulation': 'dudarev',
    }
