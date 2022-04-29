# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the :class:`aiida_quantumespresso.calculations.namelists.NamelistsCalculation` class."""
from pathlib import Path

from aiida import orm
from aiida.common import datastructures
import pytest

from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


@pytest.fixture
def generate_inputs(fixture_localhost, fixture_code, generate_remote_data):
    """Return only those inputs that the parser will expect to be there."""

    def _factory(parameters=None, settings=None, filepath_parent_folder=None):
        from aiida_quantumespresso.utils.resources import get_default_options

        inputs = {
            'code': fixture_code('quantumespresso.namelists'),
            'parameters': orm.Dict(dict=parameters or {}),
            'settings': orm.Dict(dict=settings or {}),
            'metadata': {
                'options': get_default_options()
            }
        }

        if filepath_parent_folder is not None:
            inputs['parent_folder'] = generate_remote_data(
                fixture_localhost, str(filepath_parent_folder), 'quantumespresso.pw'
            )

        return inputs

    return _factory


def test_namelists(fixture_sandbox, generate_calc_job, generate_inputs, file_regression):
    """Test a default ``NamelistCalculation``."""
    entry_point_name = 'quantumespresso.namelists'
    inputs = generate_inputs()
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    retrieve_list = ['aiida.out']
    local_copy_list = []

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.local_copy_list) == sorted(local_copy_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)
    assert calc_info.retrieve_temporary_list == []

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')


def test_namelists_parent_folder(tmp_path, fixture_sandbox, generate_calc_job, generate_inputs):
    """Test a ``NamelistCalculation`` with a parent folder that has nested directories."""
    dirpath = tmp_path / 'nested' / 'folder'
    dirpath.mkdir(parents=True)

    filepath = dirpath / 'some_file.txt'
    filepath.write_text('some context')

    entry_point_name = 'quantumespresso.namelists'
    inputs = generate_inputs(filepath_parent_folder=tmp_path)
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    cls = NamelistsCalculation
    remote = inputs['parent_folder']
    remote_path = Path(remote.get_remote_path())

    remote_copy_list = [
        (remote.computer.uuid, str(remote_path / cls._default_parent_output_folder), cls._OUTPUT_SUBFOLDER),  # pylint: disable=protected-access
    ]
    assert sorted(calc_info.remote_copy_list) == remote_copy_list
