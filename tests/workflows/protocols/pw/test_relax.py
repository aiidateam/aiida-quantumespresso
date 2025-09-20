# -*- coding: utf-8 -*-
"""Tests for the ``PwRelaxWorkChain.get_builder_from_protocol`` method."""
# pylint: disable=no-member
from aiida.engine import ProcessBuilder
import pytest

from aiida_quantumespresso.common.types import ElectronicType, RelaxType, SpinType
from aiida_quantumespresso.workflows.pw.relax import PwRelaxWorkChain

pytestmark = pytest.mark.usefixtures('pseudo_family')


def test_get_available_protocols():
    """Test ``PwRelaxWorkChain.get_available_protocols``."""
    protocols = PwRelaxWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['balanced', 'fast', 'stringent']
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``PwRelaxWorkChain.get_default_protocol``."""
    assert PwRelaxWorkChain.get_default_protocol() == 'balanced'


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

    for namespace in [builder.base_init_relax, builder.base_relax]:
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['occupations'] == 'fixed'
        assert 'degauss' not in parameters['SYSTEM']
        assert 'smearing' not in parameters['SYSTEM']


def test_spin_type(fixture_code, generate_structure):
    """Test ``PwRelaxWorkChain.get_builder_from_protocol`` with ``spin_type`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, spin_type=SpinType.COLLINEAR)

    for namespace in [builder.base_init_relax, builder.base_relax]:
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['nspin'] == 2
        assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.1}

    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, spin_type=SpinType.SPIN_ORBIT)

    for namespace in [builder.base_init_relax, builder.base_relax]:
        parameters = namespace['pw']['parameters'].get_dict()
        assert parameters['SYSTEM']['noncolin'] is True
        assert parameters['SYSTEM']['lspinorb'] is True
        assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.1}


def test_relax_type(fixture_code, generate_structure):
    """Test ``PwRelaxWorkChain.get_builder_from_protocol`` with ``spin_type`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, relax_type=RelaxType.NONE)
    assert builder.base_relax['pw']['parameters']['CONTROL']['calculation'] == 'scf'
    assert 'CELL' not in builder.base_relax['pw']['parameters'].get_dict()

    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, relax_type=RelaxType.POSITIONS)
    assert builder.base_relax['pw']['parameters']['CONTROL']['calculation'] == 'relax'
    assert 'CELL' not in builder.base_relax['pw']['parameters'].get_dict()

    with pytest.raises(ValueError):
        builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, relax_type=RelaxType.VOLUME)
        # assert builder.base['pw']['parameters']['CONTROL']['calculation'] == 'vc-relax'
        # assert builder.base['pw']['parameters']['CELL']['cell_dofree'] == 'volume'
        # assert builder.base['pw']['settings'].get_dict() == {'FIXED_COORDS': [[True, True, True], [True, True, True]]}

    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, relax_type=RelaxType.SHAPE)
    assert builder.base_relax['pw']['parameters']['CONTROL']['calculation'] == 'vc-relax'
    assert builder.base_relax['pw']['parameters']['CELL']['cell_dofree'] == 'shape'
    assert builder.base_relax['pw']['settings'].get_dict() == {'FIXED_COORDS': [[True, True, True], [True, True, True]]}

    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, relax_type=RelaxType.CELL)
    assert builder.base_relax['pw']['parameters']['CONTROL']['calculation'] == 'vc-relax'
    assert builder.base_relax['pw']['parameters']['CELL']['cell_dofree'] == 'all'
    assert builder.base_relax['pw']['settings'].get_dict() == {'FIXED_COORDS': [[True, True, True], [True, True, True]]}

    with pytest.raises(ValueError):
        builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, relax_type=RelaxType.POSITIONS_VOLUME)
        # assert builder.base['pw']['parameters']['CONTROL']['calculation'] == 'vc-relax'
        # assert builder.base['pw']['parameters']['CELL']['cell_dofree'] == 'volume'
        # assert 'settings' not in builder.base['pw']

    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, relax_type=RelaxType.POSITIONS_SHAPE)
    assert builder.base_relax['pw']['parameters']['CONTROL']['calculation'] == 'vc-relax'
    assert builder.base_relax['pw']['parameters']['CELL']['cell_dofree'] == 'shape'
    assert 'settings' not in builder.base_relax['pw']

    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, relax_type=RelaxType.POSITIONS_CELL)
    assert builder.base_relax['pw']['parameters']['CONTROL']['calculation'] == 'vc-relax'
    assert builder.base_relax['pw']['parameters']['CELL']['cell_dofree'] == 'all'
    assert 'settings' not in builder.base_relax['pw']


def test_options(fixture_code, generate_structure):
    """Test specifying ``options`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    queue_name = 'super-fast'
    withmpi = False  # The protocol default is ``True``

    options = {'queue_name': queue_name, 'withmpi': withmpi}
    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure, options=options)

    for subspace in (
        builder.base_init_relax.pw.metadata,
        builder.base_relax.pw.metadata,
    ):
        assert subspace['options']['queue_name'] == queue_name


@pytest.mark.parametrize(
    'struc_name,cell_dofree', (
        ('silicon', 'all'),
        ('2D-xy-arsenic', '2Dxy'),
        ('1D-x-carbon', 'x'),
        ('1D-y-carbon', 'y'),
        ('1D-z-carbon', 'z'),
    )
)
def test_pbc_cell(fixture_code, generate_structure, struc_name, cell_dofree):
    """Test structures with various ``pbc`` set the correct ``CELL`` parameters."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure(struc_name)

    builder = PwRelaxWorkChain.get_builder_from_protocol(code, structure)
    assert builder.base_relax.pw.parameters['CELL'].get('cell_dofree', None) == cell_dofree


@pytest.mark.parametrize(
    'overrides,warning',
    (
        # CORRECT overrides for top-level process input
        ({
            'clean_workdir': True
        }, None),
        # CORRECT overrides for nested process input
        ({
            'base_relax': {
                'kpoints_force_parity': True
            }
        }, None),
        # CORRECT overrides for nested protocol input
        ({
            'base_relax': {
                'pseudo_family': 'SSSP/1.3/PBEsol/efficiency'
            }
        }, None),
        # WRONG overrides with typo
        ({
            'clean_wokdir': True
        }, UserWarning),
        # WRONG overrides with process input at incorrect level
        ({
            'base_relax': {
                'clean_workdir': True
            }
        }, UserWarning),
        # WRONG overrides with protocol input at incorrect level
        ({
            'pseudo_family': 'SSSP/1.3/PBEsol/efficiency'
        }, UserWarning),
    )
)
def test_overrides_key_check(fixture_code, generate_structure, overrides, warning):
    """Test that the `get_builder_from_protocol()` method warns for erroneous keys in the `overrides`."""

    with pytest.warns(warning):
        PwRelaxWorkChain.get_builder_from_protocol(
            fixture_code('quantumespresso.pw'),
            generate_structure('silicon'),
            overrides=overrides,
        )
