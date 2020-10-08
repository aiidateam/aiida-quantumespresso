# -*- coding: utf-8 -*-
"""Tests for the `PwCalculation` class."""

import pytest

from aiida import orm
from aiida.common import datastructures
from aiida_quantumespresso.utils.resources import get_default_options
from aiida_quantumespresso.calculations.helpers import QEInputValidationError


def test_pw_default(fixture_sandbox, generate_calc_job, generate_inputs_pw, file_regression):
    """Test a default `PwCalculation`."""
    entry_point_name = 'quantumespresso.pw'

    inputs = generate_inputs_pw()
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    upf = inputs['pseudos']['Si']

    cmdline_params = ['-in', 'aiida.in']
    local_copy_list = [(upf.uuid, upf.filename, './pseudo/Si.upf')]
    retrieve_list = ['aiida.out', './out/aiida.save/data-file-schema.xml', './out/aiida.save/data-file.xml']
    retrieve_temporary_list = [['./out/aiida.save/K*[0-9]/eigenval*.xml', '.', 2]]

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.cmdline_params) == sorted(cmdline_params)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)
    assert sorted(calc_info.retrieve_temporary_list) == sorted(retrieve_temporary_list)
    assert sorted(calc_info.remote_symlink_list) == sorted([])

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in', 'pseudo', 'out'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')


def test_pw_ibrav(
    fixture_sandbox, generate_calc_job, fixture_code, generate_kpoints_mesh, generate_upf_data, file_regression
):
    """Test a `PwCalculation` where `ibrav` is explicitly specified."""
    entry_point_name = 'quantumespresso.pw'

    parameters = {'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutrho': 240.0, 'ecutwfc': 30.0, 'ibrav': 2}}

    # The structure needs to be rotated in the same way QE does it for ibrav=2.
    param = 5.43
    cell = [[-param / 2., 0, param / 2.], [0, param / 2., param / 2.], [-param / 2., param / 2., 0]]
    structure = orm.StructureData(cell=cell)
    structure.append_atom(position=(0., 0., 0.), symbols='Si', name='Si')
    structure.append_atom(position=(param / 4., param / 4., param / 4.), symbols='Si', name='Si')

    upf = generate_upf_data('Si')
    inputs = {
        'code': fixture_code(entry_point_name),
        'structure': structure,
        'kpoints': generate_kpoints_mesh(2),
        'parameters': orm.Dict(dict=parameters),
        'pseudos': {
            'Si': upf
        },
        'metadata': {
            'options': get_default_options()
        }
    }

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

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

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in', 'pseudo', 'out'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')


def test_pw_wrong_ibrav(fixture_sandbox, generate_calc_job, fixture_code, generate_kpoints_mesh, generate_upf_data):
    """Test that a `PwCalculation` with an incorrect `ibrav` raises."""
    entry_point_name = 'quantumespresso.pw'

    parameters = {'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutrho': 240.0, 'ecutwfc': 30.0, 'ibrav': 2}}

    # Here we use the wrong order of unit cell vectors on purpose.
    param = 5.43
    cell = [[0, param / 2., param / 2.], [-param / 2., 0, param / 2.], [-param / 2., param / 2., 0]]
    structure = orm.StructureData(cell=cell)
    structure.append_atom(position=(0., 0., 0.), symbols='Si', name='Si')
    structure.append_atom(position=(param / 4., param / 4., param / 4.), symbols='Si', name='Si')

    upf = generate_upf_data('Si')
    inputs = {
        'code': fixture_code(entry_point_name),
        'structure': structure,
        'kpoints': generate_kpoints_mesh(2),
        'parameters': orm.Dict(dict=parameters),
        'pseudos': {
            'Si': upf
        },
        'metadata': {
            'options': get_default_options()
        }
    }

    with pytest.raises(QEInputValidationError):
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)


def test_pw_ibrav_tol(fixture_sandbox, generate_calc_job, fixture_code, generate_kpoints_mesh, generate_upf_data):
    """Test that `IBRAV_TOLERANCE` controls the tolerance when checking cell consistency."""
    entry_point_name = 'quantumespresso.pw'

    parameters = {'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutrho': 240.0, 'ecutwfc': 30.0, 'ibrav': 2}}

    # The structure needs to be rotated in the same way QE does it for ibrav=2.
    param = 5.43
    eps = 0.1
    cell = [[-param / 2., eps, param / 2.], [-eps, param / 2. + eps, param / 2.], [-param / 2., param / 2., 0]]
    structure = orm.StructureData(cell=cell)
    structure.append_atom(position=(0., 0., 0.), symbols='Si', name='Si')
    structure.append_atom(position=(param / 4., param / 4., param / 4.), symbols='Si', name='Si')

    upf = generate_upf_data('Si')
    inputs = {
        'code': fixture_code(entry_point_name),
        'structure': structure,
        'kpoints': generate_kpoints_mesh(2),
        'parameters': orm.Dict(dict=parameters),
        'pseudos': {
            'Si': upf
        },
        'metadata': {
            'options': get_default_options()
        },
    }
    # Without adjusting the tolerance, the check fails.
    with pytest.raises(QEInputValidationError):
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    # After adjusting the tolerance, the input validation no longer fails.
    inputs['settings'] = orm.Dict(dict={'ibrav_cell_tolerance': eps})
    generate_calc_job(fixture_sandbox, entry_point_name, inputs)
