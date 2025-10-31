"""Tests for the `NebCalculation` class."""

import pytest
from aiida import orm

from aiida_quantumespresso.calculations.neb import NebCalculation


def test_neb_calcjob_validation(fixture_sandbox, generate_calc_job, generate_inputs):
    """Test that the `NebCalculation` still uses the `CalcJob` top-level input validation."""
    entry_point_name = 'quantumespresso.neb'

    inputs = generate_inputs(entry_point_name)
    inputs['metadata'] = {'options': {'resources': {'non-existent': 1}}}
    with pytest.raises(ValueError, match='input `metadata.options.resources` is not valid for the'):
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)


def test_parameters_validation():
    """Test the validation of the `parameters` input."""

    builder = NebCalculation.get_builder()

    parameters = {'path': {'num_of_images': 8}}

    with pytest.warns(UserWarning, match="'path' should be UPPERCASE"):
        builder.parameters = parameters

    assert builder.parameters.get_dict() == {'PATH': {'num_of_images': 8}}

    with pytest.raises(ValueError, match="'path' should be UPPERCASE"):
        builder.parameters = orm.Dict(parameters).store()
