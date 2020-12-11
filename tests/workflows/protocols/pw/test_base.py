# -*- coding: utf-8 -*-
"""Tests for the ``PwBaseWorkChain.get_builder_from_protocol`` method."""
import pytest

from aiida.engine import ProcessBuilder

from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain


def test_get_available_protocols():
    """Test ``PwBaseWorkChain.get_available_protocols``."""
    protocols = PwBaseWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['fast', 'moderate', 'precise']
    all('description' in protocol for protocol in protocols)


def test_get_default_protocol():
    """Test ``PwBaseWorkChain.get_default_protocol``."""
    assert PwBaseWorkChain.get_default_protocol() == 'moderate'


def test_default(fixture_code, generate_structure, data_regression, serialize_builder):
    """Test ``PwBaseWorkChain.get_builder_from_protocol`` for the default protocol."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()
    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure)

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))


def test_electronic_type(fixture_code, generate_structure):
    """Test ``PwBaseWorkChain.get_builder_from_protocol`` with ``electronic_type`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    with pytest.raises(NotImplementedError):
        for electronic_type in [ElectronicType.AUTOMATIC]:
            PwBaseWorkChain.get_builder_from_protocol(code, structure, electronic_type=electronic_type)

    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, electronic_type=ElectronicType.INSULATOR)
    parameters = builder.pw.parameters.get_dict()

    assert parameters['SYSTEM']['occupations'] == 'fixed'
    assert 'degauss' not in parameters['SYSTEM']
    assert 'smearing' not in parameters['SYSTEM']


def test_spin_type(fixture_code, generate_structure):
    """Test ``PwBaseWorkChain.get_builder_from_protocol`` with ``spin_type`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    with pytest.raises(NotImplementedError):
        for spin_type in [SpinType.NON_COLLINEAR, SpinType.SPIN_ORBIT]:
            PwBaseWorkChain.get_builder_from_protocol(code, structure, spin_type=spin_type)

    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, spin_type=SpinType.COLLINEAR)
    parameters = builder.pw.parameters.get_dict()

    assert parameters['SYSTEM']['nspin'] == 2
    assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.1}


@pytest.mark.parametrize('initial_magnetic_moments', ({}, {'Si1': 1.0, 'Si2': 2.0}))
def test_initial_magnetic_moments_invalid(fixture_code, generate_structure, initial_magnetic_moments):
    """Test ``PwBaseWorkChain.get_builder_from_protocol`` with invalid ``initial_magnetic_moments`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    with pytest.raises(
        ValueError, match=r'`initial_magnetic_moments` is specified but spin type `.*` is incompatible.'
    ):
        PwBaseWorkChain.get_builder_from_protocol(code, structure, initial_magnetic_moments=initial_magnetic_moments)

    with pytest.raises(ValueError):
        PwBaseWorkChain.get_builder_from_protocol(
            code, structure, initial_magnetic_moments=initial_magnetic_moments, spin_type=SpinType.COLLINEAR
        )


def test_initial_magnetic_moments(fixture_code, generate_structure):
    """Test ``PwBaseWorkChain.get_builder_from_protocol`` with ``initial_magnetic_moments`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    initial_magnetic_moments = {'Si': 1.0}
    builder = PwBaseWorkChain.get_builder_from_protocol(
        code, structure, initial_magnetic_moments=initial_magnetic_moments, spin_type=SpinType.COLLINEAR
    )
    parameters = builder.pw.parameters.get_dict()

    assert parameters['SYSTEM']['nspin'] == 2
    assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.25}
