# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the ``PwBaseWorkChain.get_builder_from_protocol`` method."""
from aiida.engine import ProcessBuilder
import pytest

from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain


def test_get_available_protocols():
    """Test ``PwBaseWorkChain.get_available_protocols``."""
    protocols = PwBaseWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['balanced', 'fast', 'stringent']
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``PwBaseWorkChain.get_default_protocol``."""
    assert PwBaseWorkChain.get_default_protocol() == 'balanced'


def test_default(fixture_code, generate_structure, data_regression, serialize_builder):
    """Test ``PwBaseWorkChain.get_builder_from_protocol`` for the default protocol."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')
    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure)

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))


@pytest.mark.parametrize('old_protocol,new_protocol', (
    ('moderate', 'balanced'),
    ('precise', 'stringent'),
))
def test_old_protocol_names(fixture_code, generate_structure, serialize_builder, old_protocol, new_protocol):
    """Test that the old protocol names still work and produce the same builder contents."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    old_builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, protocol=old_protocol)
    new_builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, protocol=new_protocol)

    assert serialize_builder(old_builder) == serialize_builder(new_builder)


def test_electronic_type(fixture_code, generate_structure):
    """Test ``PwBaseWorkChain.get_builder_from_protocol`` with ``electronic_type`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    with pytest.raises(NotImplementedError):
        for electronic_type in [ElectronicType.AUTOMATIC]:
            PwBaseWorkChain.get_builder_from_protocol(code, structure, electronic_type=electronic_type)

    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, electronic_type=ElectronicType.INSULATOR)
    parameters = builder.pw.parameters.get_dict()  # pylint: disable=no-member

    assert parameters['SYSTEM']['occupations'] == 'fixed'
    assert 'degauss' not in parameters['SYSTEM']
    assert 'smearing' not in parameters['SYSTEM']


def test_spin_type(fixture_code, generate_structure):
    """Test ``PwBaseWorkChain.get_builder_from_protocol`` with ``spin_type`` keyword."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    # Test specifying no magnetic inputs
    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure)
    assert 'starting_magnetization' not in builder.pw.parameters['SYSTEM']  # pylint: disable=no-member
    assert 'nspin' not in builder.pw.parameters['SYSTEM']  # pylint: disable=no-member

    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, spin_type=SpinType.COLLINEAR)
    parameters = builder.pw.parameters.get_dict()  # pylint: disable=no-member

    assert parameters['SYSTEM']['nspin'] == 2
    assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.1}

    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, spin_type=SpinType.SPIN_ORBIT)
    parameters = builder.pw.parameters.get_dict()  # pylint: disable=no-member

    assert parameters['SYSTEM']['noncolin'] is True
    assert parameters['SYSTEM']['lspinorb'] is True
    assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.1}


@pytest.mark.parametrize(
    'struc_name,assume_isolated', (
        ('silicon', None),
        ('2D-xy-arsenic', '2D'),
        ('1D-x-carbon', None),
        ('1D-y-carbon', None),
        ('1D-z-carbon', None),
    )
)
def test_pbc_assume_isolated(fixture_code, generate_structure, struc_name, assume_isolated):
    """Test structures with various ``pbc`` set the correct ``assume_isolated``."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure(struc_name)

    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure)
    assert builder.pw.parameters['SYSTEM'].get('assume_isolated', None) == assume_isolated


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
    parameters = builder.pw.parameters.get_dict()  # pylint: disable=no-member
    assert parameters['SYSTEM']['nspin'] == 2
    assert parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.25}


def test_overrides_pseudo_family(fixture_code, generate_structure):
    """Test specifying ``pseudo_family`` ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    overrides = {'pseudo_family': 'PseudoDojo/0.4/PBEsol/FR/standard/upf'}

    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, overrides=overrides)
    parameters = builder.pw.parameters.get_dict()
    assert parameters['SYSTEM']['ecutwfc'] == 60.0
    assert parameters['SYSTEM']['ecutrho'] == 400.0


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
    assert builder.pw.parameters['SYSTEM']['starting_magnetization'] == initial_starting_magnetization  # pylint: disable=no-member
    assert builder.pw.parameters['SYSTEM']['nspin'] == 2  # pylint: disable=no-member

    # Test that specifying `overrides` override the `initial_magnetic_moments`
    builder = PwBaseWorkChain.get_builder_from_protocol(
        code,
        structure,
        overrides=overrides,
        spin_type=SpinType.COLLINEAR,
        initial_magnetic_moments=initial_magnetic_moments
    )
    assert builder.pw.parameters['SYSTEM']['starting_magnetization'] == {'Si': 0.5}  # pylint: disable=no-member
    assert builder.pw.parameters['SYSTEM']['nspin'] == 2  # pylint: disable=no-member


def test_parameter_overrides(fixture_code, generate_structure):
    """Test specifying parameter ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    overrides = {'pw': {'parameters': {'SYSTEM': {'nbnd': 123}}}}
    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, overrides=overrides)
    assert builder.pw.parameters['SYSTEM']['nbnd'] == 123  # pylint: disable=no-member


def test_settings_overrides(fixture_code, generate_structure):
    """Test specifying settings ``overrides`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('silicon')

    overrides = {'pw': {'settings': {'cmdline': ['--kickass-mode']}}}
    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, overrides=overrides)
    assert builder.pw.settings['cmdline'] == ['--kickass-mode']  # pylint: disable=no-member


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
    metadata = builder.pw.metadata  # pylint: disable=no-member

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
    parallelization = builder.pw.parallelization  # pylint: disable=no-member

    assert parallelization['npool'] == 4
    assert parallelization['ndiag'] == 12


def test_pseudos_overrides(fixture_code, generate_structure, generate_upf_data):
    """Test specifying ``pw.pseudos`` ``overrides`` for the ``get_builder_from_protocol()`` method."""
    structure = generate_structure('cobalt-prim')
    upf_data = generate_upf_data('Co', z_valence=5.0)

    for overrides in (
        {
            'pw': {
                'pseudos': {
                    'Co': upf_data,
                },
                'parameters': {
                    'SYSTEM': {
                        'ecutrho': 123,
                        'ecutwfc': 456,
                    }
                }
            }
        },
        {
            'pseudo_family': None,
            'pw': {
                'pseudos': {
                    'Co': upf_data,
                },
                'parameters': {
                    'SYSTEM': {
                        'ecutrho': 123,
                        'ecutwfc': 456,
                    }
                }
            }
        },
    ):
        builder = PwBaseWorkChain.get_builder_from_protocol(
            code=fixture_code('quantumespresso.pw'),
            structure=structure,
            overrides=overrides,
            spin_type=SpinType.COLLINEAR
        )
        pseudos = builder.pw.pseudos  # pylint: disable=no-member
        parameters = builder.pw.parameters.get_dict()

        assert pseudos['Co'] == upf_data
        assert parameters['SYSTEM']['ecutrho'] == 123
        assert parameters['SYSTEM']['ecutwfc'] == 456
        assert parameters['SYSTEM']['starting_magnetization']['Co'] == 1.0


def test_pseudos_overrides_raises(fixture_code, generate_structure, generate_upf_data):
    """Test specifying ``pw.pseudos`` ``overrides`` for the ``get_builder_from_protocol()`` method.

    This test checks that the `get_builder_from_protocol` raises when the user overrides the `pseudos` but does not
    override the cutoffs.
    """
    structure = generate_structure('cobalt-prim')
    upf_data = generate_upf_data('Co')

    for params in (
        ({
            'pw': {
                'pseudos': {
                    'Co': upf_data
                }
            }
        }, ValueError, 'both `ecutwfc` and `ecutrho` cutoffs should be'),
        ({
            'pw': {
                'pseudos': {
                    'Si': upf_data
                }
            }
        }, ValueError, '`pseudos` override needs one value'),
        ({
            'pseudo_family': None,
            'pw': {
                'pseudos': {
                    'Co': upf_data
                }
            }
        }, ValueError, 'both `ecutwfc` and `ecutrho` cutoffs should be'),
        ({
            'pseudo_family': None,
            'pw': {
                'pseudos': {
                    'Si': upf_data
                }
            }
        }, ValueError, '`pseudos` override needs one value'),
    ):
        overrides, error_type, error_msg = params

        with pytest.raises(error_type, match=error_msg):
            PwBaseWorkChain.get_builder_from_protocol(
                code=fixture_code('quantumespresso.pw'), structure=structure, overrides=overrides
            )


def test_pseudos_family_structure_fail(fixture_code, generate_structure):
    """Test that using structure with uranium fails for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure('uranium')

    with pytest.raises(ValueError, match=r'failed to obtain recommended cutoffs for pseudo family'):
        PwBaseWorkChain.get_builder_from_protocol(
            code,
            structure,
        )


def test_options(fixture_code, generate_structure):
    """Test specifying ``options`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    queue_name = 'super-fast'
    withmpi = False  # The protocol default is ``True``

    options = {'queue_name': queue_name, 'withmpi': withmpi}
    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, options=options)
    metadata = builder.pw.metadata  # pylint: disable=no-member

    assert metadata['options']['queue_name'] == queue_name
    assert metadata['options']['withmpi'] == withmpi


def test_pop_none_overrides(fixture_code, generate_structure):
    """Test popping `None` input overrides specified in ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.pw')
    structure = generate_structure()

    overrides = {'kpoints_distance': None}
    builder = PwBaseWorkChain.get_builder_from_protocol(code, structure, overrides=overrides)

    assert 'kpoints_distance' not in builder  # pylint: disable=no-member
