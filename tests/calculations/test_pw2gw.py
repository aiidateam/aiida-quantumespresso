# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the :class:`aiida_quantumespresso.calculations.pw2gw.Pw2gwCalculation` class."""
import pathlib

from aiida import orm
from aiida.common import datastructures
import pytest

from aiida_quantumespresso.utils.resources import get_default_options


@pytest.fixture
def generate_inputs(tmp_path, fixture_localhost, fixture_code, generate_remote_data):
    """Fixture: inputs for `Pw2gwCalculation`."""

    def _factory(with_symlink=False, parameters=None):

        if parameters is None:
            parameters = {
                'INPUTPP': {
                    'qplda': False,
                    'vxcdiag': False,
                    'vkb': False,
                    'Emin': 0.0,
                    'Emax': 15.0,
                    'DeltaE': 0.001,
                }
            }

        inputs = {
            'code': fixture_code('quantumespresso.pw2gw'),
            'parent_folder': generate_remote_data(fixture_localhost, str(tmp_path), 'quantumespresso.pw'),
            'parameters': orm.Dict(parameters),
            'settings': orm.Dict({'PARENT_FOLDER_SYMLINK': with_symlink}),
            'metadata': {
                'options': get_default_options()
            }
        }
        return inputs

    return _factory


@pytest.mark.parametrize('with_symlink', (False, True))
def test_pw2gw_default(fixture_sandbox, generate_calc_job, generate_inputs, file_regression, with_symlink):
    """Test a default `Pw2gwCalculation`."""
    entry_point_name = 'quantumespresso.pw2gw'

    inputs = generate_inputs(with_symlink=with_symlink)
    remote = inputs['parent_folder']
    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    retrieve_list = ['aiida.out', 'epsX.dat', 'epsY.dat', 'epsZ.dat', 'epsTOT.dat']
    remote_copy_list = [(remote.computer.uuid, str(pathlib.Path(remote.get_remote_path()) / './pseudo/'), './pseudo/')]
    remote_symlink_list = []

    ptr = remote_copy_list if not with_symlink else remote_symlink_list
    ptr.extend([(remote.computer.uuid, str(pathlib.Path(remote.get_remote_path()) / './out/'), './out/')])

    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)
    assert sorted(calc_info.remote_copy_list) == sorted(remote_copy_list)
    assert sorted(calc_info.remote_symlink_list) == sorted(remote_symlink_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
