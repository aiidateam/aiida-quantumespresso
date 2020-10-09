# -*- coding: utf-8 -*-
"""Tests for the `PhCalculation` class."""
from aiida import orm
from aiida.common import datastructures
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.utils.resources import get_default_options

PwCalculation = CalculationFactory('quantumespresso.pw')
PhCalculation = CalculationFactory('quantumespresso.ph')


def test_ph_default(
    fixture_localhost, fixture_sandbox, generate_calc_job, fixture_code, generate_kpoints_mesh, generate_remote_data,
    file_regression
):
    """Test a default `PhCalculation`."""
    entry_point_name = 'quantumespresso.ph'
    parent_entry_point = 'quantumespresso.pw'
    remote_path = fixture_sandbox.abspath

    inputs = {
        'code': fixture_code(entry_point_name),
        'parent_folder': generate_remote_data(fixture_localhost, remote_path, parent_entry_point),
        'qpoints': generate_kpoints_mesh(2),
        'parameters': orm.Dict(dict={'INPUTPH': {}}),
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


def test_ph_qpoint_list(
    fixture_localhost, fixture_sandbox, generate_calc_job, fixture_code, generate_structure, generate_kpoints_mesh,
    generate_remote_data, file_regression
):
    """Test a `PhCalculation` with a qpoint list instead of a mesh."""
    entry_point_name = 'quantumespresso.ph'
    parent_entry_point = 'quantumespresso.pw'
    remote_path = fixture_sandbox.abspath

    structure = generate_structure()
    kpoints = generate_kpoints_mesh(2).get_kpoints_mesh(print_list=True)
    qpoints = orm.KpointsData()
    qpoints.set_cell(structure.cell)
    qpoints.set_kpoints(kpoints)

    inputs = {
        'code': fixture_code(entry_point_name),
        'parent_folder': generate_remote_data(fixture_localhost, remote_path, parent_entry_point),
        'qpoints': qpoints,
        'parameters': orm.Dict(dict={'INPUTPH': {}}),
        'metadata': {
            'options': get_default_options()
        }
    }

    generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    file_regression.check(input_written, encoding='utf-8', extension='.in')
