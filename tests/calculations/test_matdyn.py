"""Tests for the `MatdynCalculation` class."""

from pathlib import Path

from aiida import orm
from aiida.common import datastructures

from aiida_quantumespresso.calculations.matdyn import MatdynCalculation
from aiida_quantumespresso.calculations.ph import PhCalculation


def test_matdyn_default(fixture_sandbox, generate_calc_job, generate_inputs_matdyn, file_regression):
    """Test a default ``MatdynCalculation``."""
    entry_point_name = 'quantumespresso.matdyn'

    inputs = generate_inputs_matdyn()
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    force_constants = inputs['force_constants']

    local_copy_list = [(force_constants.uuid, force_constants.filename, force_constants.filename)]
    retrieve_list = ['aiida.out'] + MatdynCalculation._internal_retrieve_list

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')


def test_matdyn_elph(fixture_sandbox, generate_calc_job, generate_inputs_matdyn, file_regression):
    """Test a ``MatdynCalculation`` for el-ph properties."""
    entry_point_name = 'quantumespresso.matdyn'

    inputs = generate_inputs_matdyn(parent_folder=True)
    inputs['parameters'] = orm.Dict({'INPUT': {'la2f': True, 'dos': True}})
    kpoints = orm.KpointsData()
    kpoints.set_kpoints_mesh([2, 2, 2])
    inputs['kpoints'] = kpoints
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    force_constants = inputs['force_constants']
    local_copy_list = [(force_constants.uuid, force_constants.filename, force_constants.filename)]

    retrieve_list = ['aiida.out'] + MatdynCalculation._internal_retrieve_list
    retrieve_list += [f'a2F.dos{i}' for i in range(1, 11)]
    retrieve_list.append('lambda')
    parent_folder = inputs['parent_folder']
    dirpath = Path(parent_folder.get_remote_path()) / PhCalculation._FOLDER_ELECTRON_PHONON
    remote_copy_list = [
        (
            parent_folder.computer.uuid,
            dirpath.as_posix(),
            PhCalculation._FOLDER_ELECTRON_PHONON,
        )
    ]

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.remote_copy_list) == sorted(remote_copy_list)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
