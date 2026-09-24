"""Tests for the ``PwBandsWorkChain.get_builder_from_protocol`` method."""

import pytest
from aiida.engine import ProcessBuilder

from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain

pytestmark = pytest.mark.usefixtures('pseudo_family')


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
        PwBandsWorkChain.get_builder_from_protocol(code, structure, electronic_type=ElectronicType.AUTOMATIC)

    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure, electronic_type=ElectronicType.INSULATOR)

    for namespace in [builder.scf, builder.bands]:
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['occupations'] == 'fixed'
        assert 'degauss' not in parameters['SYSTEM']
        assert 'smearing' not in parameters['SYSTEM']


def test_spin_type(fixture_code, generate_structure):
    """Test ``PwBandsWorkChain.get_builder_from_protocol`` with ``spin_type`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure, spin_type=SpinType.COLLINEAR)

    for namespace in [builder.scf, builder.bands]:
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['nspin'] == 2
        assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.1}

    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure, spin_type=SpinType.SPIN_ORBIT)

    for namespace in [builder.scf, builder.bands]:
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['noncolin'] is True
        assert parameters['SYSTEM']['lspinorb'] is True
        assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.1}


def test_bands_kpoints_overrides(fixture_code, generate_structure, generate_kpoints_mesh):
    """Test specifying bands kpoints ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    bands_kpoints = generate_kpoints_mesh(3)
    overrides = {'bands_kpoints': bands_kpoints}
    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure, overrides=overrides)
    assert builder.bands_kpoints == bands_kpoints
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
        builder.scf.pw.metadata,
        builder.bands.pw.metadata,
    ):
        assert subspace['options']['queue_name'] == queue_name, subspace


def test_overrides_merged(fixture_code, generate_structure):
    """Test that ``overrides`` are merged onto the builder, also for dict-valued ports.

    A dict-valued input port such as ``scf.handler_overrides`` is not a namespace, so it has to be serialised to a
    ``Dict`` rather than being recursed into. This pins that these overrides land instead of being dropped or raising.
    """
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')
    overrides = {
        'nbands_factor': 5.0,
        'clean_workdir': True,
        'scf': {'handler_overrides': {'handle_out_of_walltime': {'enabled': False}}},
        'bands': {'pw': {'settings': {'cmdline': ['-nk', '4']}}},
    }
    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure, overrides=overrides)

    assert builder.nbands_factor == 5.0
    assert builder.clean_workdir
    assert builder.scf.handler_overrides.get_dict() == {'handle_out_of_walltime': {'enabled': False}}
    assert builder.bands.pw.settings.get_dict() == {'cmdline': ['-nk', '4']}
