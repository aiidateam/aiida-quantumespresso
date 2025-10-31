"""Tests for the `ProjwfcCalculation` class."""

import pytest
from aiida import orm

from aiida_quantumespresso.calculations.projwfc import ProjwfcCalculation


def test_parameters_validation():
    """Test the validation of the `parameters` input."""

    builder = ProjwfcCalculation.get_builder()

    parameters = {'projwfc': {'ngauss': 0}}

    with pytest.warns(UserWarning, match="'projwfc' should be UPPERCASE"):
        builder.parameters = parameters

    assert builder.parameters.get_dict() == {'PROJWFC': {'ngauss': 0}}

    with pytest.raises(ValueError, match="'projwfc' should be UPPERCASE"):
        builder.parameters = orm.Dict(parameters).store()
