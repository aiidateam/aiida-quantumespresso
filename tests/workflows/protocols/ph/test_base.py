# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the ``PhBaseWorkChain.get_builder_from_protocol`` method."""
from aiida.engine import ProcessBuilder
import pytest

from aiida_quantumespresso.common.types import ElectronicType
from aiida_quantumespresso.workflows.ph.base import PhBaseWorkChain


def test_get_available_protocols():
    """Test ``PhBaseWorkChain.get_available_protocols``."""
    protocols = PhBaseWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['balanced', 'fast', 'stringent']
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``PhBaseWorkChain.get_default_protocol``."""
    assert PhBaseWorkChain.get_default_protocol() == 'balanced'


def test_default(fixture_code, data_regression, serialize_builder):
    """Test ``PhBaseWorkChain.get_builder_from_protocol`` for the default protocol."""
    code = fixture_code('quantumespresso.ph')
    builder = PhBaseWorkChain.get_builder_from_protocol(code)

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))


def test_qpoints_list(fixture_code, data_regression, serialize_builder):
    """Test ``PhBaseWorkChain.get_builder_from_protocol`` qpoints in overrides for the default protocol."""
    code = fixture_code('quantumespresso.ph')
    overrides = {'qpoints': [2, 2, 2]}
    builder = PhBaseWorkChain.get_builder_from_protocol(code, overrides=overrides)
    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))


def test_electronic_type(fixture_code):
    """Test ``PhBaseWorkChain.get_builder_from_protocol`` with ``electronic_type`` keyword."""
    code = fixture_code('quantumespresso.ph')

    with pytest.raises(NotImplementedError):
        for electronic_type in [ElectronicType.AUTOMATIC]:
            PhBaseWorkChain.get_builder_from_protocol(code, electronic_type=electronic_type)

    builder = PhBaseWorkChain.get_builder_from_protocol(code, electronic_type=ElectronicType.INSULATOR)
    parameters = builder.ph.parameters.get_dict()  # pylint: disable=no-member

    assert parameters['INPUTPH']['epsil']


def test_parameter_overrides(fixture_code):
    """Test specifying parameter ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.ph')

    overrides = {'ph': {'parameters': {'INPUTHP': {'nmix_ph': 20}}}}
    builder = PhBaseWorkChain.get_builder_from_protocol(code, overrides=overrides)
    assert builder.ph.parameters['INPUTHP']['nmix_ph'] == 20  # pylint: disable=no-member


def test_settings_overrides(fixture_code):
    """Test specifying settings ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.ph')

    overrides = {'ph': {'settings': {'cmdline': ['--kickass-mode']}}}
    builder = PhBaseWorkChain.get_builder_from_protocol(code, overrides=overrides)
    assert builder.ph.settings['cmdline'] == ['--kickass-mode']  # pylint: disable=no-member


def test_metadata_overrides(fixture_code):
    """Test specifying metadata ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.ph')

    overrides = {'ph': {'metadata': {'options': {'resources': {'num_machines': 1e90}, 'max_wallclock_seconds': 1}}}}
    builder = PhBaseWorkChain.get_builder_from_protocol(code, overrides=overrides)
    metadata = builder.ph.metadata  # pylint: disable=no-member

    assert metadata['options']['resources']['num_machines'] == 1e90
    assert metadata['options']['max_wallclock_seconds'] == 1


def test_options(fixture_code):
    """Test specifying ``options`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.ph')

    queue_name = 'super-fast'
    withmpi = False  # The protocol default is ``True``

    options = {'queue_name': queue_name, 'withmpi': withmpi}
    builder = PhBaseWorkChain.get_builder_from_protocol(code, options=options)
    metadata = builder.ph.metadata  # pylint: disable=no-member

    assert metadata['options']['queue_name'] == queue_name
    assert metadata['options']['withmpi'] == withmpi


def test_parent_folder(fixture_code, generate_remote_data, fixture_localhost, fixture_sandbox):
    """Test specifying ``options`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.ph')
    remote_folder = generate_remote_data(fixture_localhost, fixture_sandbox.abspath, 'quantumespresso.pw')

    builder = PhBaseWorkChain.get_builder_from_protocol(code, parent_folder=remote_folder)

    assert builder.ph.parent_folder == remote_folder  # pylint: disable=no-member
