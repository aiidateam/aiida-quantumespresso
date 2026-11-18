"""Tests for the `NebCalculation` class."""

import pytest
from aiida import orm

from aiida_quantumespresso.calculations.neb import NebCalculation


def test_neb_calcjob_validation(fixture_sandbox, generate_calc_job, generate_inputs):
    """Test that the `NebCalculation` still uses the `CalcJob` top-level input validation."""
    entry_point_name = 'quantumespresso.neb'

    inputs = generate_inputs(entry_point_name)
    inputs['metadata'] = {'options': {'resources': {'non-existent': 1}}}
    with pytest.raises(ValueError, match='input `metadata.options.resources` is not valid for the'):
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)


def test_parameters_validation():
    """Test the validation of the `parameters` input."""

    builder = NebCalculation.get_builder()

    parameters = {'path': {'num_of_images': 8}}

    with pytest.warns(UserWarning, match="'path' should be UPPERCASE"):
        builder.parameters = parameters

    assert builder.parameters.get_dict() == {'PATH': {'num_of_images': 8}}

    with pytest.raises(ValueError, match="'path' should be UPPERCASE"):
        builder.parameters = orm.Dict(parameters).store()


@pytest.fixture
def generate_trajectory_custom_kinds():
    """Return a ``TrajectoryData`` of two images whose two iron sites have the custom kind names ``Fe1`` and ``Fe2``."""

    def _generate_trajectory_custom_kinds():
        cell = [[2.87, 0, 0], [0, 2.87, 0], [0, 0, 2.87]]
        structures = []
        for shift in (0.0, 0.1):
            structure = orm.StructureData(cell=cell)
            structure.append_atom(position=(shift, 0.0, 0.0), symbols='Fe', name='Fe1')
            structure.append_atom(position=(1.435, 1.435, 1.435), symbols='Fe', name='Fe2')
            structures.append(structure)
        return orm.TrajectoryData(structures)

    return _generate_trajectory_custom_kinds


def test_neb_custom_kind_names(
    fixture_sandbox, generate_calc_job, generate_inputs, generate_upf_data, generate_trajectory_custom_kinds
):
    """Test that kind names that are not chemical symbols survive the ``images`` round trip.

    This is e.g. the case for collinear antiferromagnets, where the magnetic sublattices are distinct kinds of the same
    element, and ``starting_magnetization`` is assigned per kind.
    """
    entry_point_name = 'quantumespresso.neb'

    inputs = generate_inputs(entry_point_name)
    inputs['images'] = generate_trajectory_custom_kinds()
    pseudo = generate_upf_data('Fe')
    inputs['pw']['pseudos'] = {'Fe1': pseudo, 'Fe2': pseudo}
    inputs['pw']['parameters'] = orm.Dict(
        {
            'CONTROL': {'calculation': 'scf'},
            'SYSTEM': {
                'ecutrho': 240.0,
                'ecutwfc': 30.0,
                'nspin': 2,
                'starting_magnetization': {'Fe1': 0.5, 'Fe2': -0.5},
            },
        }
    )

    generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    for filename in ('pw_1.in', 'pw_2.in'):
        with fixture_sandbox.open(filename) as handle:
            content = handle.read()
        assert 'ntyp = 2' in content
        assert 'starting_magnetization(1) =   5.0000000000d-01' in content
        assert 'starting_magnetization(2) =  -5.0000000000d-01' in content

        species = content.split('ATOMIC_SPECIES')[1].split('ATOMIC_POSITIONS')[0].split()
        assert species[0::3] == ['Fe1', 'Fe2']
        assert species[2::3] == ['Fe.upf', 'Fe.upf']

        positions = content.split('ATOMIC_POSITIONS')[1].split('K_POINTS')[0].splitlines()[1:]
        assert [line.split()[0] for line in positions if line.strip()] == ['Fe1', 'Fe2']


def test_neb_pseudos_kinds_mismatch(fixture_sandbox, generate_calc_job, generate_inputs, generate_upf_data):
    """Test that a mismatch between the ``pseudos`` and the kinds of the ``images`` is reported as such."""
    entry_point_name = 'quantumespresso.neb'

    inputs = generate_inputs(entry_point_name)
    inputs['pw']['pseudos'] = {'He': generate_upf_data('He')}

    with pytest.raises(ValueError, match='Mismatch between the defined pseudos and the list of kinds'):
        generate_calc_job(fixture_sandbox, entry_point_name, inputs)
