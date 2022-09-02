# -*- coding: utf-8 -*-
"""Tests for the `Q2rCalculation` class."""
# pylint: disable=protected-access
from pathlib import Path

from aiida.common import datastructures

from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.calculations.q2r import Q2rCalculation


def test_q2r_default(fixture_sandbox, generate_calc_job, generate_inputs_q2r, file_regression):
    """Test a default `Q2rCalculation`."""
    entry_point_name = 'quantumespresso.q2r'

    inputs = generate_inputs_q2r()
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    parent_folder = inputs['parent_folder']
    remote_copy_folder = Path(parent_folder.get_remote_path()) / PhCalculation._FOLDER_DYNAMICAL_MATRIX

    remote_copy_list = [(parent_folder.computer.uuid, str(remote_copy_folder), PhCalculation._FOLDER_DYNAMICAL_MATRIX)]
    retrieve_list = [Q2rCalculation._DEFAULT_OUTPUT_FILE] + Q2rCalculation._internal_retrieve_list

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.remote_copy_list) == sorted(remote_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
