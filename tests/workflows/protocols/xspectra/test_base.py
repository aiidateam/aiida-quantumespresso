# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the ``XspectraBaseWorkChain.get_builder_from_protocol`` method."""
from aiida.engine import ProcessBuilder
from aiida.plugins import WorkflowFactory
import pytest

XspectraBaseWorkChain = WorkflowFactory('quantumespresso.xspectra.base')


@pytest.fixture
def get_xspectra_generator_inputs(fixture_code, generate_structure):
    """Generate a set of default inputs for the ``XspectraBaseWorkChain.get_builder_from_protocol()`` method."""
    return {
        'pw_code': fixture_code('quantumespresso.pw'),
        'xs_code': fixture_code('quantumespresso.xspectra'),
        'structure': generate_structure('silicon'),
    }


def test_get_available_protocols():
    """Test ``XspectraBaseWorkChain.get_available_protocols()``."""
    protocols = XspectraBaseWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['full', 'half', 'none', 'xch_fix', 'xch_smear']
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``XspectraBaseWorkChain.get_default_protocol``."""
    assert XspectraBaseWorkChain.get_default_protocol() == 'full'


def test_default(get_xspectra_generator_inputs, data_regression, serialize_builder):
    """Test ``XspectraBaseWorkChain.get_builder_from_protocol`` for the default protocol."""
    inputs = get_xspectra_generator_inputs
    builder = XspectraBaseWorkChain.get_builder_from_protocol(**inputs)

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))
