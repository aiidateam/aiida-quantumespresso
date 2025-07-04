# -*- coding: utf-8 -*-
"""Tests for the ``PwBandsWorkChain.get_builder_from_protocol`` method."""
from aiida.engine import ProcessBuilder
import pytest

from aiida_quantumespresso.common.types import ElectronicType, RelaxType, SpinType
from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain


def test_get_available_protocols():
    """Test ``PwBandsWorkChain.get_available_protocols``."""
    protocols = PwBandsWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['balanced', 'fast', 'stringent']
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``PwBandsWorkChain.get_default_protocol``."""
    assert PwBandsWorkChain.get_default_protocol() == 'balanced'


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


def test_bands_kpoints_overrides(fixture_code, generate_structure, generate_kpoints_mesh):
    """Test specifying bands kpoints ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    bands_kpoints = generate_kpoints_mesh(3)
    overrides = {'bands_kpoints': bands_kpoints}
    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure, overrides=overrides)
    assert builder.bands_kpoints == bands_kpoints  # pylint: disable=no-member
    assert 'bands_kpoints_distance' not in builder


def test_options(fixture_code, generate_structure):
    """Test specifying ``options`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    queue_name = 'super-fast'
    withmpi = False  # The protocol default is ``True``

    options = {'queue_name': queue_name, 'withmpi': withmpi}
    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure, options=options)

    for subspace in (
        builder.relax.base.pw.metadata,
        builder.scf.pw.metadata,  # pylint: disable=no-member
        builder.bands.pw.metadata,  # pylint: disable=no-member
    ):
        assert subspace['options']['queue_name'] == queue_name, subspace
