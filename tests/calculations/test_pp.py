# -*- coding: utf-8 -*-
"""Tests for the `PpCalculation` class."""
from __future__ import absolute_import

from aiida import orm
from aiida.common import datastructures
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.utils.resources import get_default_options

PwCalculation = CalculationFactory('quantumespresso.pw')
PpCalculation = CalculationFactory('quantumespresso.pp')


def test_pp_default(
    aiida_profile, fixture_localhost, fixture_sandbox, generate_calc_job, fixture_code, generate_remote_data,
    file_regression
):
    """Test a default `PpCalculation`."""
    entry_point_name = 'quantumespresso.pp'
    parent_entry_point = 'quantumespresso.pw'
    remote_path = fixture_sandbox.abspath

    inputs = {
        'code': fixture_code(entry_point_name),
        'parent_folder': generate_remote_data(fixture_localhost, remote_path, parent_entry_point),
        'parameters': orm.Dict(dict={
            'INPUTPP': {
                'plot_num': 1
            },
            'PLOT': {
                'iflag': 3
            }
        }),
        'metadata': {
            'options': get_default_options()
        }
    }

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    retrieve_list = ['aiida.out']
    retrieve_temporary_list = ['aiida.fileout']
    local_copy_list = []

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)
    assert sorted(calc_info.retrieve_temporary_list) == sorted(retrieve_temporary_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')


def test_pp_keep_plot_file(
    aiida_profile, fixture_localhost, fixture_sandbox, generate_calc_job, fixture_code, generate_remote_data,
    file_regression
):
    """Test a `PpCalculation` where we want to retrieve the plot file."""
    entry_point_name = 'quantumespresso.pp'
    parent_entry_point = 'quantumespresso.pw'
    remote_path = fixture_sandbox.abspath

    inputs = {
        'code': fixture_code(entry_point_name),
        'parent_folder': generate_remote_data(fixture_localhost, remote_path, parent_entry_point),
        'parameters': orm.Dict(dict={
            'INPUTPP': {
                'plot_num': 1
            },
            'PLOT': {
                'iflag': 3
            }
        }),
        'metadata': {
            'options': {
                'keep_plot_file': True
            },
        }
    }

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    retrieve_list = ['aiida.out', 'aiida.fileout']
    retrieve_temporary_list = []
    local_copy_list = []

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)
    assert sorted(calc_info.retrieve_temporary_list) == sorted(retrieve_temporary_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
