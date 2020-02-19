# -*- coding: utf-8 -*-
"""Tests for the `EpwCalculation` class."""
from __future__ import absolute_import

from aiida import orm
from aiida.common import datastructures
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.utils.resources import get_default_options

PwCalculation = CalculationFactory('quantumespresso.pw')
PhCalculation = CalculationFactory('quantumespresso.ph')
EpwCalculation = CalculationFactory('quantumespresso.epw')

def test_example():
    assert 1==1

def test_epw_default(
    aiida_profile, fixture_localhost, fixture_sandbox, generate_calc_job, fixture_code, generate_kpoints_mesh,
    generate_remote_data, file_regression
):
    """Test a default `EpwCalculation`."""
    entry_point_name = 'quantumespresso.epw'
    parent_entry_point = 'quantumespresso.pw'
    parent_entry_point_ph = 'quantumespresso.ph'
    remote_path = fixture_sandbox.abspath

    inputs = {
        'code': fixture_code(entry_point_name),
        'parent_folder': generate_remote_data(fixture_localhost, remote_path, parent_entry_point),
        'parent_folder_ph': generate_remote_data(fixture_localhost, remote_path, parent_entry_point_ph),
        'qpoints': generate_kpoints_mesh(2),
        'kpoints': generate_kpoints_mesh(6),
        'parameters': orm.Dict(dict={'INPUTEPW': {}}),
        'metadata': {
            'options': get_default_options()
        }
    }

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    cmdline_params = ['-in', 'aiida.in']
    retrieve_list = ['./out/_ph0/aiida.phsave/tensors.xml', 'DYN_MAT', 'aiida.out']
    local_copy_list = []

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.codes_info[0].cmdline_params) == sorted(cmdline_params)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)
    assert sorted(calc_info.remote_symlink_list) == sorted([])

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['DYN_MAT', 'aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')


