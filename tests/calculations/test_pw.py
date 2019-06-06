# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `PwCalculation` class."""
from __future__ import absolute_import

from aiida import orm
from aiida.common import datastructures

from aiida_quantumespresso.utils.resources import get_default_options


def test_pw_default(fixture_database, fixture_computer_localhost, fixture_sandbox_folder, generate_calc_job,
    generate_code_localhost, generate_structure, generate_kpoints_mesh, generate_upf_data, file_regression):
    """Test a default `PwCalculation`."""
    entry_point_name = 'quantumespresso.pw'

    parameters = {
        'CONTROL': {
            'calculation': 'scf'
        },
        'SYSTEM': {
            'ecutrho': 240.0,
            'ecutwfc': 30.0
        }
    }

    upf = generate_upf_data('Si')
    inputs = {
        'code': generate_code_localhost(entry_point_name, fixture_computer_localhost),
        'structure': generate_structure('Si'),
        'kpoints': generate_kpoints_mesh(2),
        'parameters': orm.Dict(dict=parameters),
        'pseudos': {'Si': upf},
        'metadata': {'options': get_default_options()}
    }

    calc_info = generate_calc_job(fixture_sandbox_folder, entry_point_name, inputs)

    cmdline_params = ['-in', 'aiida.in']
    local_copy_list = [(upf.uuid, upf.filename, u'./pseudo/Si.upf')]
    retrieve_list = ['aiida.out', './out/aiida.save/data-file-schema.xml', './out/aiida.save/data-file.xml']
    retrieve_temporary_list = [['./out/aiida.save/K*[0-9]/eigenval*.xml', '.', 2]]

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.cmdline_params) == sorted(cmdline_params)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)
    assert sorted(calc_info.retrieve_temporary_list) == sorted(retrieve_temporary_list)
    assert sorted(calc_info.remote_symlink_list) == sorted([])

    with fixture_sandbox_folder.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox_folder.get_content_list()) == sorted(['aiida.in', 'pseudo', 'out'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
