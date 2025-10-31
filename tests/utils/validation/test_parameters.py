"""Tests for parameter casing validation."""

import pytest
from aiida import orm

from aiida_quantumespresso.utils.validation.parameters import validate_parameters


@pytest.mark.parametrize(
    'input_params,expected_params,should_warn',
    [
        # Correct casing - no warning
        (
            {'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutwfc': 30.0}},
            {'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutwfc': 30.0}},
            False,
        ),
        # Lowercase namelists - warning and correction
        (
            {'control': {'calculation': 'scf'}, 'system': {'ecutwfc': 30.0}},
            {'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutwfc': 30.0}},
            True,
        ),
        # Uppercase parameters - warning and correction
        (
            {'CONTROL': {'CALCULATION': 'scf'}, 'SYSTEM': {'ECUTWFC': 30.0}},
            {'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutwfc': 30.0}},
            True,
        ),
        # Mixed case - warning and correction
        (
            {'control': {'CALCULATION': 'scf'}, 'SYSTEM': {'ECUTWFC': 30.0}},
            {'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutwfc': 30.0}},
            True,
        ),
    ],
)
def test_validate_unstored_node(input_params, expected_params, should_warn):
    """Test validation on unstored nodes with various casing."""
    params = orm.Dict(input_params)

    if should_warn:
        with pytest.warns(UserWarning):
            result = validate_parameters(params, None)
    else:
        result = validate_parameters(params, None)

    assert result is None
    assert params.get_dict() == expected_params


def test_validate_stored_node_incorrect_casing_returns_error():
    """Test that incorrect casing on stored nodes returns an error message."""
    params = orm.Dict({'control': {'calculation': 'scf'}})
    params.store()

    result = validate_parameters(params, None)

    assert result is not None
    assert "'control' should be UPPERCASE" in result


def test_validate_stored_node_correct_casing_passes():
    """Test that correct casing on stored nodes passes."""
    params = orm.Dict({'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutwfc': 30.0}})
    params.store()

    assert validate_parameters(params, None) is None


def test_validate_namelist_value_not_dict():
    """Test that non-dictionary namelist values return an error."""
    params = orm.Dict({'CONTROL': 'not a dictionary'})

    result = validate_parameters(params, None)

    assert result is not None
    assert "'CONTROL' should contain a dictionary of parameters" in result
