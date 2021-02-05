# -*- coding: utf-8 -*-
"""Tests for the ``PwBandsWorkChain.get_builder_from_protocol`` method."""
import pytest

from aiida.engine import ProcessBuilder

from aiida_quantumespresso.common.types import ElectronicType, RelaxType, SpinType
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


def test_electronic_type(fixture_code, generate_structure):
    """Test ``PwBandsWorkChain.get_builder_from_protocol`` with ``electronic_type`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    with pytest.raises(NotImplementedError):
        for electronic_type in [ElectronicType.AUTOMATIC]:
            PwBandsWorkChain.get_builder_from_protocol(code, structure, electronic_type=electronic_type)

    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure, electronic_type=ElectronicType.INSULATOR)

    for namespace in [builder.relax['base'], builder.scf, builder.bands]:
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['occupations'] == 'fixed'
        assert 'degauss' not in parameters['SYSTEM']
        assert 'smearing' not in parameters['SYSTEM']


def test_spin_type(fixture_code, generate_structure):
    """Test ``PwBandsWorkChain.get_builder_from_protocol`` with ``spin_type`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    with pytest.raises(NotImplementedError):
        for spin_type in [SpinType.NON_COLLINEAR, SpinType.SPIN_ORBIT]:
            PwBandsWorkChain.get_builder_from_protocol(code, structure, spin_type=spin_type)

    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure, spin_type=SpinType.COLLINEAR)

    for namespace in [builder.relax['base'], builder.scf, builder.bands]:
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['nspin'] == 2
        assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.1}


def test_relax_type(fixture_code, generate_structure):
    """Test ``PwBandsWorkChain.get_builder_from_protocol`` setting the ``relax_type`` input."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure, relax_type=RelaxType.NONE)
    assert builder.relax['base']['pw']['parameters']['CONTROL']['calculation'] == 'scf'
    assert 'CELL' not in builder.relax['base']['pw']['parameters'].get_dict()
