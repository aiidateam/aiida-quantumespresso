# -*- coding: utf-8 -*-
"""Tests for the ``PwRelaxWorkChain.get_builder_from_protocol`` method."""
import pytest

from aiida.engine import ProcessBuilder

from aiida_quantumespresso.common.types import ElectronicType
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


def test_electronic_type(fixture_code, generate_structure):
    """Test ``PwRelaxWorkChain.get_builder_from_protocol`` with ``electronic_type`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    with pytest.raises(NotImplementedError):
        for electronic_type in [ElectronicType.AUTOMATIC]:
            PwRelaxWorkChain.get_builder_from_protocol(code, structure, electronic_type=electronic_type)

    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, electronic_type=ElectronicType.INSULATOR)

    for namespace in [builder.base, builder.base_final_scf]:
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['occupations'] == 'fixed'
        assert 'degauss' not in parameters['SYSTEM']
        assert 'smearing' not in parameters['SYSTEM']
