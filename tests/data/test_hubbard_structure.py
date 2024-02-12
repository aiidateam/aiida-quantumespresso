# -*- coding: utf-8 -*-
"""Tests for the :mod:`data.hubbard_structure` module."""
# pylint: disable=redefined-outer-name,protected-access
import pytest

from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData


@pytest.fixture
def generate_hubbard():
    """Return a `Hubbard` instance."""

    def _generate_hubbard():
        from aiida_quantumespresso.common.hubbard import Hubbard
        return Hubbard.from_list([(0, '1s', 0, '1s', 5.0, (0, 0, 0), 'Ueff')])

    return _generate_hubbard


@pytest.fixture
def generate_hubbard_structure(generate_structure):
    """Return a `HubbardStructureData` instance."""

    def _generate_hubbard_structure():
        from aiida_quantumespresso.common.hubbard import Hubbard
        structure = generate_structure('silicon-kinds')
        hp_list = [(0, '1s', 0, '1s', 5.0, (0, 0, 0), 'Ueff')]
        hubbard = Hubbard.from_list(hp_list)
        return HubbardStructureData.from_structure(structure=structure, hubbard=hubbard)

    return _generate_hubbard_structure


@pytest.mark.usefixtures('aiida_profile')
def test_valid_init(generate_hubbard):
    """Test the constructor."""
    cell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    sites = [['Si', 'Si0', (0, 0, 0)]]
    hubbard = generate_hubbard()
    hubbard_structure = HubbardStructureData(cell=cell, sites=sites, hubbard=hubbard)
    assert hubbard_structure.cell == cell
    assert hubbard_structure.kinds[0].symbol == sites[0][0]
    assert hubbard_structure.sites[0].kind_name == sites[0][1]
    assert hubbard_structure.sites[0].position == sites[0][2]
    assert hubbard_structure.hubbard == hubbard


@pytest.mark.usefixtures('aiida_profile')
def test_from_structure(generate_structure, generate_hubbard):
    """Test the `from_structure` method."""
    structure = generate_structure()
    hubbard = generate_hubbard()
    hubbard_structure = HubbardStructureData.from_structure(structure=structure, hubbard=hubbard)
    assert hubbard_structure.cell == structure.cell
    assert len(hubbard_structure.sites) == len(structure.sites)
    assert hubbard_structure.hubbard == hubbard  # sanity check

    structure = generate_structure('silicon-kinds')
    hubbard_structure = HubbardStructureData.from_structure(structure=structure)
    assert hubbard_structure.get_site_kindnames() == structure.get_site_kindnames()
    assert len(hubbard_structure.kinds) == 2


@pytest.mark.usefixtures('aiida_profile')
@pytest.mark.parametrize('structure_name', ('silicon',))
@pytest.mark.parametrize(
    'parameters', (
        ((0, '1s', 0, '1s', 5.0, (0, 0, 0), 'Ueff'),),
        ((0, '1s', 1, '1s', 5.0, None, 'V'),),
        (
            (0, '1s', 0, '1s', 5.0, (0, 0, 0), 'Ueff'),
            (0, '1s', 0, '1s', 5.0, (0, 0, 0), 'Ueff'),
        ),
        (
            (0, '1s', 0, '1s', 5.0, (0, 0, 0), 'Ueff'),
            (0, '1s', 1, '1s', 5.0, (0, 1, 0), 'V'),
        ),
    )
)
def test_append_hubbard_parameters(data_regression, generate_structure, structure_name, parameters):
    """Test the `append_hubbard_parameters` method."""
    hubbard_structure = HubbardStructureData.from_structure(generate_structure(structure_name))

    for parameter in parameters:
        hubbard_structure.append_hubbard_parameter(*parameter)

    data_regression.check(hubbard_structure.hubbard.to_list())
    assert len(hubbard_structure.hubbard.parameters) == len(set(parameters))


@pytest.mark.parametrize('parameter', (
    (0, '1s', 1, '1s', 5.0, None, 'V'),
    (0, '1s', 1, '1s', 5.0, (0, 0, 0), 'V'),
))
def test_append_hubbard_parameters_invalid_index(generate_structure, parameter):
    """Test the `append_hubbard_parameters` method with invalid index."""
    hubbard_structure = HubbardStructureData.from_structure(generate_structure('cobalt-prim'))

    with pytest.raises(ValueError, match='atom_index and neighbour_index must be within the range'):
        hubbard_structure.append_hubbard_parameter(*parameter)


@pytest.mark.usefixtures('aiida_profile')
def test_pop_hubbard_parameters(generate_hubbard_structure):
    """Test the `pop_hubbard_parameters` method."""
    hubbard_structure = generate_hubbard_structure()
    hubbard_structure.pop_hubbard_parameters(0)
    assert len(hubbard_structure.hubbard.parameters) == 0


@pytest.mark.usefixtures('aiida_profile')
def test_clear_hubbard_parameters(generate_hubbard_structure):
    """Test the `clear_hubbard_parameters` method."""
    hubbard_structure = generate_hubbard_structure()
    hubbard_structure.clear_hubbard_parameters()
    assert len(hubbard_structure.hubbard.parameters) == 0


@pytest.mark.usefixtures('aiida_profile')
def test_is_storable(generate_hubbard_structure):
    """Test the storing does not throw errors."""
    hubbard_structure = generate_hubbard_structure()
    hubbard_structure.store()
    assert hubbard_structure.is_stored


@pytest.mark.usefixtures('aiida_profile')
def test_initialize_intersites_hubbard(generate_hubbard_structure):
    """Test the `initialize_intersites_hubbard` method."""
    hubbard_structure = generate_hubbard_structure()
    hubbard_structure.initialize_intersites_hubbard('Si', '1s', 'Si', '2s', 0, 'Ueff', False)

    # !WARNING! This is not the expected behavior, as we would like it to initialize
    #           intersites among first neighbours. The method was designed for different
    #           interacting species. We may want to improve it for this special cases,
    #           although it is still rather futuristic.
    assert hubbard_structure.hubbard.parameters[1].to_tuple() == (0, '1s', 0, '2s', 0.0, (0, 0, 0), 'Ueff')

    hubbard_structure.clear_hubbard_parameters()
    hubbard_structure.initialize_intersites_hubbard('Si0', '1s', 'Si1', '2s', 0.0, 'Ueff')
    assert (0, '1s', 1, '2s', 0, (-1, 0, 0), 'Ueff') in hubbard_structure.hubbard.to_list()
    assert len(hubbard_structure.hubbard.parameters) == 1

    with pytest.raises(ValueError):
        hubbard_structure.initialize_intersites_hubbard('Mg', '1s', 'Si1', '2s', 0, 'Ueff')


@pytest.mark.usefixtures('aiida_profile')
def test_initialize_onsites_hubbard(generate_hubbard_structure):
    """Test the `initialize_onsites_hubbard` method."""
    hubbard_structure = generate_hubbard_structure()

    hubbard_structure.clear_hubbard_parameters()
    hubbard_structure.initialize_onsites_hubbard('Si', '1s', 0.0, 'Ueff', False)

    assert (0, '1s', 0, '1s', 0, (0, 0, 0), 'Ueff') in hubbard_structure.hubbard.to_list()
    assert len(hubbard_structure.hubbard.parameters) == 2

    hubbard_structure.clear_hubbard_parameters()
    hubbard_structure.initialize_onsites_hubbard('Si0', '1s', 0.0, 'Ueff', True)

    assert len(hubbard_structure.hubbard.parameters) == 1


@pytest.mark.usefixtures('aiida_profile')
def test_get_one_kind_index(generate_hubbard_structure):
    """Test the `_get_one_kind_index` method."""
    hubbard_structure = generate_hubbard_structure()
    assert hubbard_structure._get_one_kind_index('Si0') == [0]
    assert hubbard_structure._get_one_kind_index('Si1') == [1]


@pytest.mark.usefixtures('aiida_profile')
def test_get_symbol_indices(generate_hubbard_structure):
    """Test the `_get_symbol_indices` method."""
    hubbard_structure = generate_hubbard_structure()
    assert hubbard_structure._get_symbol_indices('Si') == [0, 1]
