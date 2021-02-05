# -*- coding: utf-8 -*-
"""Tests for the `CpCalculation` class."""
import pytest

from aiida.common import datastructures


@pytest.mark.parametrize('autopilot', [True, False])
def test_cp_autopilot(fixture_sandbox, generate_calc_job, generate_inputs_cp, file_regression, autopilot):
    """Test a default `CpCalculation`."""
    entry_point_name = 'quantumespresso.cp'

    inputs = generate_inputs_cp(autopilot=autopilot)
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    upf = inputs['pseudos']['Si']

    cmdline_params = ['-in', 'aiida.in']
    local_copy_list = [(upf.uuid, upf.filename, './pseudo/Si.upf')]
    retrieve_list = [
        'aiida.out', './out/aiida_51.save/data-file-schema.xml', './out/aiida_51.save/data-file.xml', './out/aiida.cel',
        './out/aiida.con', './out/aiida.eig', './out/aiida.evp', './out/aiida.for', './out/aiida.nos',
        './out/aiida.pol', './out/aiida.pos', './out/aiida.spr', './out/aiida.str', './out/aiida.the',
        './out/aiida.vel', './out/aiida.wfc', './out/aiida_51.save/print_counter',
        './out/aiida_51.save/print_counter.xml'
    ]
    retrieve_temporary_list = [['./out/aiida.save/K*[0-9]/eigenval*.xml', '.', 2]]

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.codes_info[0].cmdline_params) == sorted(cmdline_params)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)
    assert sorted(calc_info.retrieve_temporary_list) == sorted(retrieve_temporary_list)
    assert sorted(calc_info.remote_symlink_list) == sorted([])

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in', 'pseudo', 'out'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
