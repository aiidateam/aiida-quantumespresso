# -*- coding: utf-8 -*-
"""Tests for the `PwCalculation` class."""

import pytest

from aiida import orm
from aiida.common import datastructures
from aiida.common.warnings import AiidaDeprecationWarning
from aiida.common.exceptions import InputValidationError
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
    assert isinstance(calc_info.codes_info[0], datastructures.CodeInfo)
    assert sorted(calc_info.codes_info[0].cmdline_params) == cmdline_params
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
    assert isinstance(calc_info.codes_info[0], datastructures.CodeInfo)
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


def test_pw_parallelization_inputs(fixture_sandbox, generate_calc_job, generate_inputs_pw):
    """Test that the parallelization settings are set correctly in the commandline params."""
    entry_point_name = 'quantumespresso.pw'

    inputs = generate_inputs_pw()
    inputs['parallelization'] = orm.Dict(dict={'npool': 4, 'nband': 2, 'ntg': 3, 'ndiag': 12})
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    cmdline_params = ['-npool', '4', '-nband', '2', '-ntg', '3', '-ndiag', '12', '-in', 'aiida.in']

    # Check that the command-line parameters are as expected.
    assert calc_info.codes_info[0].cmdline_params == cmdline_params


@pytest.mark.parametrize('flag_name', ['npool', 'nk', 'nband', 'nb', 'ntg', 'nt', 'northo', 'ndiag', 'nd'])
def test_pw_parallelization_deprecation(fixture_sandbox, generate_calc_job, generate_inputs_pw, flag_name):
    """Test the deprecation warning on specifying parallelization flags manually.

    Test that passing parallelization flags in the `settings['CMDLINE']
    emits an `AiidaDeprecationWarning`.
    """
    entry_point_name = 'quantumespresso.pw'

    inputs = generate_inputs_pw()
    extra_cmdline_args = [f'-{flag_name}', '2']
    inputs['settings'] = orm.Dict(dict={'CMDLINE': extra_cmdline_args})
    with pytest.warns(AiidaDeprecationWarning) as captured_warnings:
        calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
        assert calc_info.codes_info[0].cmdline_params == extra_cmdline_args + ['-in', 'aiida.in']
    assert any('parallelization flags' in str(warning.message) for warning in captured_warnings.list)


def test_pw_parallelization_conflict_error(fixture_sandbox, generate_calc_job, generate_inputs_pw):
    """Test conflict between `settings['CMDLINE']` and `parallelization`.

    Test that passing the same parallelization flag (modulo aliases)
    manually in `settings['CMDLINE']` and in the `parallelization`
    input raises an `InputValidationError`.
    """
    entry_point_name = 'quantumespresso.pw'

    inputs = generate_inputs_pw()
    extra_cmdline_args = ['-nk', '2']
    inputs['settings'] = orm.Dict(dict={'CMDLINE': extra_cmdline_args})
    inputs['parallelization'] = orm.Dict(dict={'npool': 2})
    with pytest.raises(InputValidationError) as exc:
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    assert 'conflicts' in str(exc.value)


def test_pw_parallelization_incorrect_flag(fixture_sandbox, generate_calc_job, generate_inputs_pw):
    """Test that passing a non-existing parallelization flag raises.

    Test that specifying an non-existing parallelization flag in
    the `parallelization` `Dict` raises a `ValueError`.
    """
    entry_point_name = 'quantumespresso.pw'

    inputs = generate_inputs_pw()
    inputs['parallelization'] = orm.Dict(dict={'invalid_flag_name': 2})
    with pytest.raises(ValueError) as exc:
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    assert 'Unknown' in str(exc.value)


def test_pw_parallelization_incorrect_value(fixture_sandbox, generate_calc_job, generate_inputs_pw):
    """Test that passing a non-integer parallelization flag raises.

    Test that specifying an non-integer parallelization flag value in
    the `parallelization` `Dict` raises a `ValueError`.
    """
    entry_point_name = 'quantumespresso.pw'

    inputs = generate_inputs_pw()
    inputs['parallelization'] = orm.Dict(dict={'npool': 2.2})
    with pytest.raises(ValueError) as exc:
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    assert 'integer' in str(exc.value)


def test_pw_parallelization_duplicate_cmdline_flag(fixture_sandbox, generate_calc_job, generate_inputs_pw):
    """Test that passing two different aliases to the same parallelization flag raises."""
    entry_point_name = 'quantumespresso.pw'

    inputs = generate_inputs_pw()
    inputs['settings'] = orm.Dict(dict={'CMDLINE': ['-nk', '2', '-npools', '2']})
    with pytest.raises(InputValidationError) as exc:
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    assert 'Conflicting' in str(exc.value)
