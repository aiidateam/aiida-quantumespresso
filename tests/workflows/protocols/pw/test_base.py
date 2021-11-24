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
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``PwBaseWorkChain.get_default_protocol``."""
    assert PwBaseWorkChain.get_default_protocol() == 'moderate'


def test_default(fixture_code, generate_structure, data_regression, serialize_builder):
    """Test ``PwBaseWorkChain.get_builder_from_protocol`` for the default protocol."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')
    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure)

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))


def test_electronic_type(fixture_code, generate_structure):
    """Test ``PwBaseWorkChain.get_builder_from_protocol`` with ``electronic_type`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

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
    structure = generate_structure('silicon')

    # Test specifying no magnetic inputs
    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure)
    assert 'starting_magnetization' not in builder.pw.parameters['SYSTEM']
    assert 'nspin' not in builder.pw.parameters['SYSTEM']

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
    structure = generate_structure('silicon')

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
    structure = generate_structure('silicon')

    initial_magnetic_moments = {'Si': 1.0}
    builder = PwBaseWorkChain.get_builder_from_protocol(
        code, structure, initial_magnetic_moments=initial_magnetic_moments, spin_type=SpinType.COLLINEAR
    )
    parameters = builder.pw.parameters.get_dict()
    assert parameters['SYSTEM']['nspin'] == 2
    assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.25}


def test_magnetization_overrides(fixture_code, generate_structure):
    """Test magnetization ``overrides`` for the ``PwBaseWorkChain.get_builder_from_protocol`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')
    initial_magnetic_moments = {'Si': 1.0}
    initial_starting_magnetization = {'Si': 0.5}
    overrides = {'pw': {'parameters': {'SYSTEM': {'starting_magnetization': initial_starting_magnetization}}}}

    # Test specifying `starting_magnetization` via the `overrides`
    builder = PwBaseWorkChain.get_builder_from_protocol(
        code, structure, overrides=overrides, spin_type=SpinType.COLLINEAR
    )
    assert builder.pw.parameters['SYSTEM']['starting_magnetization'] == initial_starting_magnetization
    assert builder.pw.parameters['SYSTEM']['nspin'] == 2

    # Test that specifying `overrides` override the `initial_magnetic_moments`
    builder = PwBaseWorkChain.get_builder_from_protocol(
        code,
        structure,
        overrides=overrides,
        spin_type=SpinType.COLLINEAR,
        initial_magnetic_moments=initial_magnetic_moments
    )
    assert builder.pw.parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.5}
    assert builder.pw.parameters['SYSTEM']['nspin'] == 2


def test_parameter_overrides(fixture_code, generate_structure):
    """Test specifying parameter ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    overrides = {'pw': {'parameters': {'SYSTEM': {'nbnd': 123}}}}
    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, overrides=overrides)
    assert builder.pw.parameters['SYSTEM']['nbnd'] == 123


def test_settings_overrides(fixture_code, generate_structure):
    """Test specifying settings ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    overrides = {'pw': {'settings': {'cmdline': ['--kickass-mode']}}}
    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, overrides=overrides)
    assert builder.pw.settings['cmdline'] == ['--kickass-mode']


def test_metadata_overrides(fixture_code, generate_structure):
    """Test specifying metadata ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    overrides = {'pw': {'metadata': {'options': {'resources': {'num_machines': 1e90}, 'max_wallclock_seconds': 1}}}}
    builder = PwBaseWorkChain.get_builder_from_protocol(
        code,
        structure,
        overrides=overrides,
    )
    metadata = builder.pw.metadata

    assert metadata['options']['resources']['num_machines'] == 1e90
    assert metadata['options']['max_wallclock_seconds'] == 1


def test_parallelization_overrides(fixture_code, generate_structure):
    """Test specifying parallelization ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    overrides = {'pw': {'parallelization': {'npool': 4, 'ndiag': 12}}}
    builder = PwBaseWorkChain.get_builder_from_protocol(
        code,
        structure,
        overrides=overrides,
    )
    parallelization = builder.pw.parallelization

    assert parallelization['npool'] == 4
    assert parallelization['ndiag'] == 12


def test_pseudos_overrides(fixture_code, generate_structure, generate_upf_data):
    """Test specifying ``pw.pseudos`` ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')
    silicon_pseudo = generate_upf_data('Si')

    overrides = {'pw': {'pseudos': {'Si': silicon_pseudo}}}
    builder = PwBaseWorkChain.get_builder_from_protocol(
        code,
        structure,
        overrides=overrides,
    )
    pseudos = builder.pw.pseudos

    assert pseudos['Si'] == silicon_pseudo


def test_pseudos_family_structure_fail(fixture_code, generate_structure):
    """Test that using structure with uranium fails for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('uranium')

    with pytest.raises(ValueError, match=r'failed to obtain recommended cutoffs for pseudo family'):
        PwBaseWorkChain.get_builder_from_protocol(
            code,
            structure,
        )
