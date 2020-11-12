# -*- coding: utf-8 -*-
"""Tests for the ``PwRelaxWorkChain.get_builder_from_protocol`` method."""
from aiida.engine import ProcessBuilder
from aiida_quantumespresso.workflows.pw.relax import PwRelaxWorkChain


def test_get_available_protocols():
    """Test ``PwRelaxWorkChain.get_available_protocols``."""
    protocols = PwRelaxWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['fast', 'moderate', 'precise']
    all('description' in protocol for protocol in protocols)


def test_get_default_protocol():
    """Test ``PwRelaxWorkChain.get_default_protocol``."""
    assert PwRelaxWorkChain.get_default_protocol() == 'moderate'


def test_default(fixture_code, generate_structure, data_regression, serialize_builder):
    """Test ``PwRelaxWorkChain.get_builder_from_protocol`` for the default protocol."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()
    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure)

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))
