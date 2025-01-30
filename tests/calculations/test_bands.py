# -*- coding: utf-8 -*-
"""Tests for the `BandsCalculation` class."""
# pylint: disable=protected-access
from aiida.common import datastructures

from aiida_quantumespresso.calculations.bands import BandsCalculation


def test_bands_default(fixture_sandbox, generate_calc_job, generate_inputs_bands, file_regression):
    """Test a default `BandsCalculation`."""
    entry_point_name = 'quantumespresso.bands'

    inputs = generate_inputs_bands()
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    retrieve_list = [BandsCalculation._DEFAULT_OUTPUT_FILE] + BandsCalculation._internal_retrieve_list

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
