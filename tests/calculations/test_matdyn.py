# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `MatdynCalculation` class."""
from __future__ import absolute_import

import os

from aiida.common import datastructures
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.data.force_constants import ForceConstantsData
from aiida_quantumespresso.utils.resources import get_default_options

Q2rCalculation = CalculationFactory('quantumespresso.q2r')
MatdynCalculation = CalculationFactory('quantumespresso.matdyn')


def test_matdyn_default(fixture_database, fixture_computer_localhost, fixture_sandbox_folder, generate_calc_job,
    generate_code_localhost, generate_structure, generate_kpoints_mesh, file_regression):
    """Test a default `MatdynCalculation`."""
    entry_point_name = 'quantumespresso.matdyn'

    filepath = os.path.join(os.path.dirname(__file__), 'fixtures', 'matdyn', 'default', 'force_constants.dat')
    force_constants = ForceConstantsData(filepath)

    inputs = {
        'code': generate_code_localhost(entry_point_name, fixture_computer_localhost),
        'force_constants': force_constants,
        'kpoints': generate_kpoints_mesh(2),
        'metadata': {'options': get_default_options()}
    }

    calc_info = generate_calc_job(fixture_sandbox_folder, entry_point_name, inputs)

    local_copy_list = [(force_constants.uuid, force_constants.filename, force_constants.filename)]
    retrieve_list = ['aiida.out'] + MatdynCalculation._internal_retrieve_list

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox_folder.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox_folder.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
