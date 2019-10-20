# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `ProjwfcCalculation` class."""
from __future__ import absolute_import

from aiida import orm
from aiida.common import datastructures

from aiida_quantumespresso.utils.resources import get_default_options


def test_projwfc_default(
    fixture_database, fixture_computer_localhost, fixture_sandbox_folder, generate_calc_job, generate_code_localhost,
    file_regression
):
    """Test a default `ProjwfcCalculation`."""
    entry_point_name = 'quantumespresso.projwfc'

    parameters = {'PROJWFC': {'emin': -1, 'emax': 1, 'DeltaE': 0.01, 'ngauss': 0, 'degauss': 0.01}}
    parent_folder = orm.RemoteData(computer=fixture_computer_localhost, remote_path='path/on/remote')
    parent_folder.store()

    inputs = {
        'code': generate_code_localhost(entry_point_name, fixture_computer_localhost),
        'parameters': orm.Dict(dict=parameters),
        'parent_folder': parent_folder,
        'metadata': {
            'options': get_default_options()
        }
    }

    calc_info = generate_calc_job(fixture_sandbox_folder, entry_point_name, inputs)
    code_info = calc_info.codes_info[0]

    # Check the attributes of the returned `CodeInfo`
    assert isinstance(code_info, datastructures.CodeInfo)
    assert code_info.stdin_name == 'aiida.in'
    assert code_info.stdout_name == 'aiida.out'
    assert code_info.cmdline_params == []
    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.local_copy_list) == sorted([])
    assert sorted(calc_info.remote_copy_list
                  ) == sorted([(parent_folder.computer.uuid, 'path/on/remote/./out/', './out/')])
    assert sorted(calc_info.retrieve_list) == sorted(['aiida.pdos*'])
    assert sorted(calc_info.retrieve_temporary_list) == sorted(['aiida.out'])

    with fixture_sandbox_folder.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox_folder.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
