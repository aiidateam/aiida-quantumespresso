# -*- coding: utf-8 -*-
"""Tests for the `MatdynCalculation` class."""
from aiida.common import datastructures
from aiida.plugins import CalculationFactory

MatdynCalculation = CalculationFactory('quantumespresso.matdyn')


def test_matdyn_default(fixture_sandbox, generate_calc_job, generate_inputs_matdyn, file_regression):
    """Test a default `MatdynCalculation`."""
    entry_point_name = 'quantumespresso.matdyn'

    inputs = generate_inputs_matdyn()
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    force_constants = inputs['force_constants']

    local_copy_list = [(force_constants.uuid, force_constants.filename, force_constants.filename)]
    retrieve_list = ['aiida.out'] + MatdynCalculation._internal_retrieve_list  # pylint: disable=protected-access

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
