"""Tests for the `Q2rCalculation` class."""

from pathlib import Path

import pytest
from aiida import orm
from aiida.common import datastructures

from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.calculations.q2r import Q2rCalculation


def test_q2r_default(fixture_sandbox, generate_calc_job, generate_inputs_q2r, file_regression):
    """Test a default ``Q2rCalculation``."""
    entry_point_name = 'quantumespresso.q2r'

    inputs = generate_inputs_q2r()
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)
    parent_folder = inputs['parent_folder']
    remote_copy_folder = Path(parent_folder.get_remote_path()) / PhCalculation._FOLDER_DYNAMICAL_MATRIX

    remote_copy_list = [
        (
            parent_folder.computer.uuid,
            str(remote_copy_folder),
            PhCalculation._FOLDER_DYNAMICAL_MATRIX,
        )
    ]
    retrieve_list = [Q2rCalculation._DEFAULT_OUTPUT_FILE] + Q2rCalculation._internal_retrieve_list

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.remote_copy_list) == sorted(remote_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')


@pytest.mark.parametrize(
    ('key', 'value'),
    [
        ('la2F', True),
        ('zasr', 'simple'),
    ],
)
def test_q2r_parameters(fixture_sandbox, generate_calc_job, generate_inputs_q2r, file_regression, key, value):
    """Test the ``parameters`` input for the ``Q2rCalculation`` plugin."""
    entry_point_name = 'quantumespresso.q2r'

    inputs = generate_inputs_q2r()
    inputs['parameters'] = orm.Dict({'INPUT': {key: value}})
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    parent_folder = inputs['parent_folder']
    remote_copy_folder = Path(parent_folder.get_remote_path()) / PhCalculation._FOLDER_DYNAMICAL_MATRIX
    remote_copy_list = [
        (
            parent_folder.computer.uuid,
            str(remote_copy_folder),
            PhCalculation._FOLDER_DYNAMICAL_MATRIX,
        )
    ]
    retrieve_list = [Q2rCalculation._DEFAULT_OUTPUT_FILE] + Q2rCalculation._internal_retrieve_list

    # Test that `elph_dir` is added properly to the remote_copy_list
    if inputs['parameters']['INPUT'].get('la2F', None) is True:
        dirpath = Path(parent_folder.get_remote_path()) / PhCalculation._FOLDER_ELECTRON_PHONON
        remote_copy_list.append(
            (
                parent_folder.computer.uuid,
                dirpath.as_posix(),
                PhCalculation._FOLDER_ELECTRON_PHONON,
            )
        )

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.remote_copy_list) == sorted(remote_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
