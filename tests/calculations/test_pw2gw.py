# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the `Pw2gwCalculation` class."""
import os
import pytest

from aiida import orm
from aiida.common import datastructures


@pytest.fixture()
def parameters():
    """Fixture: parameters for pw calculation."""
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

    return parameters


@pytest.fixture(params=[
    False,
], ids=[
    'base',
])
def settings(request):
    """Fixture: parameters for z2pack calculation."""
    settings = {}

    if request.param:
        settings['PARENT_FOLDER_SYMLINK'] = True

    return settings


@pytest.fixture
def remote(fixture_localhost, tmpdir, generate_remote_data):
    """Fixture: Remote folder created by a CalcJob with inputs."""
    remote = generate_remote_data(
        fixture_localhost,
        str(tmpdir),
        'quantumespresso.pw',
    )

    return remote


@pytest.fixture()
def inputs(fixture_code, remote, parameters, settings):
    """Fixture: inputs for Z2packBaseWorkChain."""
    from aiida_quantumespresso.utils.resources import get_default_options

    inputs = {
        'code': fixture_code('quantumespresso.pw2gw'),
        'parent_folder': remote,
        'parameters': orm.Dict(dict=parameters),
        'settings': orm.Dict(dict=settings),
        'metadata': {
            'options': get_default_options()
        }
    }
    return inputs


@pytest.mark.parametrize(
    'settings,with_symlink', [(False, False), (True, True)], ids=['base', 'with_symlink'], indirect=['settings']
)
def test_pw2gw_default(fixture_sandbox, generate_calc_job, file_regression, remote, inputs, with_symlink):
    """Test a default `Pw2gwCalculation`."""
    entry_point_name = 'quantumespresso.pw2gw'

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    retrieve_list = ['aiida.out', 'epsX.dat', 'epsY.dat', 'epsZ.dat', 'epsTOT.dat']

    remote_copy_list = [
        (remote.computer.uuid, os.path.join(remote.get_remote_path(), path), path) for path in ['./pseudo/']
    ]
    remote_symlink_list = []

    ptr = remote_copy_list if not with_symlink else remote_symlink_list
    ptr.extend([(remote.computer.uuid, os.path.join(remote.get_remote_path(), path), path) for path in ['./out/']])

    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)
    assert sorted(calc_info.remote_copy_list) == sorted(remote_copy_list)
    assert sorted(calc_info.remote_symlink_list) == sorted(remote_symlink_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
