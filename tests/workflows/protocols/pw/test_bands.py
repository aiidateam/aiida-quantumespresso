# -*- coding: utf-8 -*-
"""Tests for the ``PwBandsWorkChain.get_builder_from_protocol`` method."""
from aiida.engine import ProcessBuilder
from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain


def test_get_available_protocols():
    """Test ``PwBandsWorkChain.get_available_protocols``."""
    protocols = PwBandsWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['fast', 'moderate', 'precise']
    all('description' in protocol for protocol in protocols)


def test_get_default_protocol():
    """Test ``PwBandsWorkChain.get_default_protocol``."""
    assert PwBandsWorkChain.get_default_protocol() == 'moderate'


def test_default(fixture_code, generate_structure, data_regression, serialize_builder):
    """Test ``PwBandsWorkChain.get_builder_from_protocol`` for the default protocol."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()
    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure)

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))
