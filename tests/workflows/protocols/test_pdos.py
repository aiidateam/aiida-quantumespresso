# -*- coding: utf-8 -*-
"""Tests for the ``PwBandsWorkChain.get_builder_from_protocol`` method."""
import pytest

from aiida.engine import ProcessBuilder
from aiida.plugins import WorkflowFactory

from aiida_quantumespresso.common.types import ElectronicType, SpinType

PdosWorkChain = WorkflowFactory('quantumespresso.pdos')


def test_get_available_protocols():
    """Test ``PdosWorkChain.get_available_protocols``."""
    protocols = PdosWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['fast', 'moderate', 'precise']
    all('description' in protocol for protocol in protocols)


def test_get_default_protocol():
    """Test ``PdosWorkChain.get_default_protocol``."""
    assert PdosWorkChain.get_default_protocol() == 'moderate'


def test_default(fixture_code, generate_structure, data_regression, serialize_builder):
    """Test ``PdosWorkChain.get_builder_from_protocol`` for the default protocol."""
    pw_code = fixture_code('quantumespresso.pw')
    dos_code = fixture_code('quantumespresso.dos')
    projwfc_code = fixture_code('quantumespresso.projwfc')
    structure = generate_structure()
    builder = PdosWorkChain.get_builder_from_protocol(
        pw_code=pw_code, dos_code=dos_code, projwfc_code=projwfc_code, structure=structure
    )

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))


def test_electronic_type(fixture_code, generate_structure):
    """Test ``PdosWorkChain.get_builder_from_protocol`` with ``electronic_type`` keyword."""
    pw_code = fixture_code('quantumespresso.pw')
    dos_code = fixture_code('quantumespresso.dos')
    projwfc_code = fixture_code('quantumespresso.projwfc')
    structure = generate_structure()

    with pytest.raises(NotImplementedError):
        builder = PdosWorkChain.get_builder_from_protocol(
            pw_code=pw_code,
            dos_code=dos_code,
            projwfc_code=projwfc_code,
            structure=structure,
            electronic_type=ElectronicType.AUTOMATIC
        )

    builder = PdosWorkChain.get_builder_from_protocol(
        pw_code=pw_code,
        dos_code=dos_code,
        projwfc_code=projwfc_code,
        structure=structure,
        electronic_type=ElectronicType.INSULATOR
    )
    for namespace, occupations in zip((builder.scf, builder.nscf), ('fixed', 'tetrahedra')):
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['occupations'] == occupations
        assert 'degauss' not in parameters['SYSTEM']
        assert 'smearing' not in parameters['SYSTEM']


def test_spin_type(fixture_code, generate_structure):
    """Test ``PdosWorkChain.get_builder_from_protocol`` with ``spin_type`` keyword."""
    pw_code = fixture_code('quantumespresso.pw')
    dos_code = fixture_code('quantumespresso.dos')
    projwfc_code = fixture_code('quantumespresso.projwfc')
    structure = generate_structure()

    with pytest.raises(NotImplementedError):
        for spin_type in [SpinType.NON_COLLINEAR, SpinType.SPIN_ORBIT]:
            builder = PdosWorkChain.get_builder_from_protocol(
                pw_code=pw_code, dos_code=dos_code, projwfc_code=projwfc_code, structure=structure, spin_type=spin_type
            )
    builder = PdosWorkChain.get_builder_from_protocol(
        pw_code=pw_code,
        dos_code=dos_code,
        projwfc_code=projwfc_code,
        structure=structure,
        spin_type=SpinType.COLLINEAR
    )
    for namespace in [builder.scf, builder.nscf]:
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['nspin'] == 2
        assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.1}
