# -*- coding: utf-8 -*-
"""Tests for :py:mod:`~aiida_quantumespresso.common.hubbard`."""
from pydantic import ValidationError
import pytest

from aiida_quantumespresso.common.hubbard import Hubbard, HubbardParameters, HubbardProjectors

# Maybe we should define a Class for testing?
valid_parameters = {
    'atom_index': 0,
    'atom_manifold': '3d',
    'neighbour_index': 1,
    'neighbour_manifold': '2p',
    'translation_vector': [0, 0, 0],
    'hubbard_value': 5.0,
    'hubbard_type': 'Dudarev-U',
}


@pytest.fixture
def get_hubbard_parameters():
    """Return an `HubbardParameters` intstance."""

    def _get_hubbard_parameters(overrides=None):
        """Return an `HubbardParameters` intstance."""
        inputs = valid_parameters.copy()

        if overrides:
            inputs.update(overrides)

        return HubbardParameters(**inputs)

    return _get_hubbard_parameters


@pytest.fixture
def get_hubbard():
    """Return an `Hubbard` intstance."""

    def _get_hubbard():
        """Return an `Hubbard` intstance."""
        hp = HubbardParameters(**valid_parameters)

        return Hubbard(hubbard_parameters=[hp, hp])

    return _get_hubbard


def test_valid_hubbard_projectors():
    """Test valid inputs for py:meth:`HubbardProjectors`."""
    from aiida_quantumespresso.common.hubbard import allowed_projectors
    for projectors in allowed_projectors:
        hubbard_projectors = HubbardProjectors(hubbard_projectors=projectors)
        assert hubbard_projectors.hubbard_projectors == projectors


def test_invalid_hubbard_projectors():
    """Test valid inputs for py:meth:`HubbardProjectors`."""
    invalid_name = 'invalid-name'
    with pytest.raises(ValidationError):
        HubbardProjectors(hubbard_projectors=invalid_name)


def test_safe_hubbard_parameters(get_hubbard_parameters):
    """Test valid inputs are stored correctly for py:meth:`HubbardParameters`."""
    hp_dict = get_hubbard_parameters().dict()
    assert hp_dict == valid_parameters


def test_from_to_list_parameters(get_hubbard_parameters):
    """Test py:meth:`HubbardParameters.to_list` and py:meth:`HubbardParameters.from_list`."""
    hp = get_hubbard_parameters()
    hp_list = [0, '3d', 1, '2p', 5.0, [0, 0, 0], 'Dudarev-U']
    assert hp.to_list() == hp_list
    hp = HubbardParameters.from_list(hp_list)
    assert hp.dict() == valid_parameters


@pytest.mark.parametrize(
    'overrides',
    [
        {
            'atom_index': 0
        },  # making sure is valid - some validators (e.g. PositiveInt) do not accept it
        {
            'atom_manifold': '3d-2p'
        },
        {
            'translation_vector': [0, -1, +1]
        },  # this type also valid
        {
            'hubbard_type': 'Liechtenstein-B'
        }
    ]
)
def test_valid_hubbard_parameters(get_hubbard_parameters, overrides):
    """Test valid inputs for py:meth:`HubbardParameters`."""
    hp_dict = get_hubbard_parameters(overrides=overrides).dict()
    new_dict = valid_parameters.copy()
    new_dict.update(overrides)
    assert hp_dict == new_dict


@pytest.mark.parametrize(
    'overrides', [
        {
            'atom_index': -1
        },
        {
            'atom_index': 0.5
        },
        {
            'atom_manifold': '3z'
        },
        {
            'atom_manifold': '3d2p'
        },
        {
            'atom_manifold': '3d-3p-2s'
        },
        {
            'translation_vector': [0, 0]
        },
        {
            'translation_vector': [0, 0, 0, 0]
        },
        {
            'translation_vector': [0, 0, -1.5]
        },
        {
            'hubbard_type': 'liechtenstein-L'
        },
        {
            'hubbard_type': 'hubbard-L'
        },
        {
            'hubbard_type': 'DudarevU'
        },
    ]
)
def test_invalid_hubbard_parameters(get_hubbard_parameters, overrides):
    """Test invalid inputs for py:meth:`HubbardParameters`."""
    with pytest.raises(ValidationError):
        get_hubbard_parameters(overrides=overrides)
    # maybe we should add the check on the Exception info as well


def test_from_to_list_hubbard(get_hubbard):
    """Test py:meth:`Hubbard.to_list` and py:meth:`Hubbard.from_list`."""
    hubbard = get_hubbard()
    hp_list = [0, '3d', 1, '2p', 5.0, [0, 0, 0], 'Dudarev-U']

    hubbard_list = [hp_list, hp_list]
    assert hubbard.to_list() == hubbard_list

    hubbard = Hubbard.from_list(hubbard_list)
    assert hubbard.dict() == {
        'hubbard_parameters': [valid_parameters, valid_parameters],
        'hubbard_projectors': 'ortho-atomic',
    }
