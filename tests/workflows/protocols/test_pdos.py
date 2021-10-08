# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the ``PdosWorkChain.get_builder_from_protocol`` method."""
import pytest

from aiida.engine import ProcessBuilder
from aiida.plugins import WorkflowFactory

from aiida_quantumespresso.common.types import ElectronicType, SpinType

PdosWorkChain = WorkflowFactory('quantumespresso.pdos')


@pytest.fixture
def get_pdos_generator_inputs(fixture_code, generate_structure):
    """Generate a set of default inputs for the ``PdosWorkChain.get_builder_from_protocol()`` method."""
    return {
        'pw_code': fixture_code('quantumespresso.pw'),
        'dos_code': fixture_code('quantumespresso.dos'),
        'projwfc_code': fixture_code('quantumespresso.projwfc'),
        'structure': generate_structure('silicon')
    }


def test_get_available_protocols():
    """Test ``PdosWorkChain.get_available_protocols``."""
    protocols = PdosWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['fast', 'moderate', 'precise']
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``PdosWorkChain.get_default_protocol``."""
    assert PdosWorkChain.get_default_protocol() == 'moderate'


def test_default(get_pdos_generator_inputs, data_regression, serialize_builder):
    """Test ``PdosWorkChain.get_builder_from_protocol`` for the default protocol."""
    inputs = get_pdos_generator_inputs
    builder = PdosWorkChain.get_builder_from_protocol(**inputs)

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))


def test_electronic_type(get_pdos_generator_inputs):
    """Test ``PdosWorkChain.get_builder_from_protocol`` with ``electronic_type`` keyword."""
    with pytest.raises(NotImplementedError):
        builder = PdosWorkChain.get_builder_from_protocol(
            **get_pdos_generator_inputs, electronic_type=ElectronicType.AUTOMATIC
        )
    builder = PdosWorkChain.get_builder_from_protocol(
        **get_pdos_generator_inputs, electronic_type=ElectronicType.INSULATOR
    )
    for namespace, occupations in zip((builder.scf, builder.nscf), ('fixed', 'tetrahedra')):
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['occupations'] == occupations
        assert 'degauss' not in parameters['SYSTEM']
        assert 'smearing' not in parameters['SYSTEM']


def test_spin_type(get_pdos_generator_inputs):
    """Test ``PdosWorkChain.get_builder_from_protocol`` with ``spin_type`` keyword."""
    with pytest.raises(NotImplementedError):
        for spin_type in [SpinType.NON_COLLINEAR, SpinType.SPIN_ORBIT]:
            builder = PdosWorkChain.get_builder_from_protocol(**get_pdos_generator_inputs, spin_type=spin_type)
    builder = PdosWorkChain.get_builder_from_protocol(**get_pdos_generator_inputs, spin_type=SpinType.COLLINEAR)
    for namespace in [builder.scf, builder.nscf]:
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['nspin'] == 2
        assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.1}


def test_nscf_smearing_raises(get_pdos_generator_inputs):
    """Test ``PdosWorkChain.get_builder_from_protocol`` fails when NSCF uses smearing."""
    overrides = {'nscf': {'pw': {'parameters': {'SYSTEM': {'occupations': 'smearing'}}}}}
    with pytest.raises(ValueError, match=r'`SYSTEM.occupations` in `nscf.pw.parameters`'):
        PdosWorkChain.get_builder_from_protocol(**get_pdos_generator_inputs, overrides=overrides)
