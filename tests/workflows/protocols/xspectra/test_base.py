# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the ``XspectraBaseWorkChain.get_builder_from_protocol`` method."""
from aiida.engine import ProcessBuilder
from aiida.plugins import WorkflowFactory
import pytest

XspectraBaseWorkChain = WorkflowFactory('quantumespresso.xspectra.base')


@pytest.fixture
def get_xspectra_generator_inputs(fixture_code, generate_structure, generate_singlefile_data):
    """Generate a set of default inputs for the ``XspectraBaseWorkChain.get_builder_from_protocol()`` method."""
    return {
        'code': fixture_code('quantumespresso.xspectra'),
        'structure': generate_structure('silicon'),
        'core_wfc_data': generate_singlefile_data(),
        'abs_atom_marker': 'Si',
    }


def test_get_available_protocols():
    """Test ``XspectraBaseWorkChain.get_available_protocols()``."""
    protocols = XspectraBaseWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['one_shot', 'production', 're_plot']
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``XspectraBaseWorkChain.get_default_protocol``."""
    assert XspectraBaseWorkChain.get_default_protocol() == 'one_shot'


def test_default(get_xspectra_generator_inputs, data_regression, serialize_builder):
    """Test ``XspectraBaseWorkChain.get_builder_from_protocol`` for the default protocol."""
    inputs = get_xspectra_generator_inputs
    builder = XspectraBaseWorkChain.get_builder_from_protocol(**inputs)

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))


def test_options(get_xspectra_generator_inputs):
    """Test the ``options`` input for the XspectraBaseWorkChain.get_builder_from_protocol() method."""
    queue_name = 'super-fast'
    withmpi = False  # The protocol default is ``True``

    options = {'queue_name': queue_name, 'withmpi': withmpi}
    builder = XspectraBaseWorkChain.get_builder_from_protocol(**get_xspectra_generator_inputs, options=options)

    for subspace in (builder.xspectra.metadata,):
        assert subspace['options']['queue_name'] == queue_name, subspace
