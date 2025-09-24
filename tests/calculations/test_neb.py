# -*- coding: utf-8 -*-
"""Tests for the `NebCalculation` class."""
import pytest


def test_neb_calcjob_validation(fixture_sandbox, generate_calc_job, generate_inputs):
    """Test that the `NebCalculation` still uses the `CalcJob` top-level input validation."""
    entry_point_name = 'quantumespresso.neb'

    inputs = generate_inputs(entry_point_name)
    inputs['metadata'] = {'options': {'resources': {'non-existent': 1}}}
    with pytest.raises(ValueError, match='input `metadata.options.resources` is not valid for the'):
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)
