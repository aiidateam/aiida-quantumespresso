# -*- coding: utf-8 -*-
"""Tests for the :mod:`utils.hubbard` module."""
# pylint: disable=redefined-outer-name
import os

from aiida.orm import StructureData
from ase.io import read
import numpy as np
import pytest

from aiida_quantumespresso.common.hubbard import Hubbard
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiida_quantumespresso.utils.hubbard import HubbardUtils, initialize_hubbard_parameters


@pytest.fixture
def generate_hubbard_structure():
    """Return a `HubbardStructureData` instance."""

    def _generate_hubbard_structure(parameters=None, projectors=None, formulation=None):
        """Return a `HubbardStructureData` instance."""
        a, b, c, d = 1.40803, 0.81293, 4.68453, 1.62585  # pylint: disable=invalid-name
        cell = [[a, -b, c], [0.0, d, c], [-a, -b, c]]  # pylint: disable=invalid-name
        positions = [[0, 0, 0], [0, 0, 3.6608], [0, 0, 10.392], [0, 0, 7.0268]]
        symbols = ['Co', 'O', 'O', 'Li']
        structure = StructureData(cell=cell)
        for position, symbol in zip(positions, symbols):
            structure.append_atom(position=position, symbols=symbol)

        if parameters is None:
            parameters = [(0, '3d', 0, '3d', 5.0, (0, 0, 0), 'U'), (0, '3d', 1, '2p', 1.0, (0, 0, 0), 'V')]
        if projectors is None:
            projectors = 'ortho-atomic'
        if formulation is None:
            formulation = 'dudarev'
        args = (parameters, projectors, formulation)
        hubbard = Hubbard.from_list(*args)

        return HubbardStructureData.from_structure(structure=structure, hubbard=hubbard)

    return _generate_hubbard_structure


@pytest.fixture
def generate_hubbard_utils(generate_hubbard_structure):
    """Return a `HubbardUtils` instance."""

    def _generate_hubbard_utils(**kwargs):
        return HubbardUtils(hubbard_structure=generate_hubbard_structure(**kwargs))

    return _generate_hubbard_utils


@pytest.fixture
def filepath_hubbard(filepath_tests):
    """Return the absolute filepath to the directory containing the file `fixtures`."""
    return os.path.join(filepath_tests, 'utils', 'fixtures', 'hubbard')


@pytest.mark.usefixtures('aiida_profile')
def test_valid_init(generate_hubbard_utils, generate_hubbard_structure):
    """Test the constructor."""
    hubbard_utils = generate_hubbard_utils()
    assert hubbard_utils.hubbard_structure.hubbard == generate_hubbard_structure().hubbard


@pytest.mark.usefixtures('aiida_profile')
def test_is_intersite(generate_hubbard_structure):
    """Test the `is_intersite_hubbard` method."""
    from aiida_quantumespresso.utils.hubbard import is_intersite_hubbard
    assert is_intersite_hubbard(generate_hubbard_structure().hubbard)


@pytest.mark.parametrize(('filename', 'projectors'), (
    ('HUBBARD.dat', 'ortho-atomic'),
    ('HUBBARD_2.dat', 'atomic'),
    ('HUBBARD_3.dat', 'atomic'),
))
@pytest.mark.usefixtures('aiida_profile')
def test_invertibility(generate_hubbard_utils, filepath_hubbard, filename, projectors):
    """Test the invertibility of the `get_hubbard_card` and `parse_hubbard_card` method.

    We test both methods at once, first parsing and then matching the generated card
    with the raw parsed file. While this may seem trivial, the parsing does not store
    the indecis, but the translation vectors.
    Thus, mathematically we are testing ~ F^-1[F(x)] = x
    """
    hubbard_utils = generate_hubbard_utils()
    filepath = os.path.join(filepath_hubbard, filename)
    hubbard_utils.parse_hubbard_dat(filepath)

    hubbard = hubbard_utils.hubbard_structure.hubbard
    assert hubbard.projectors == projectors
    assert hubbard.formulation == 'dudarev'

    hubbard_data = []
    with open(filepath, encoding='utf-8') as file:
        lines = file.readlines()
        for line in lines:
            if line.strip().split()[0] != '#':
                hubbard_data.append(tuple(line.strip().split()))

    hubbard_data.pop(0)  # removing header

    parsed_data = []
    card = hubbard_utils.get_hubbard_card()
    card = card.splitlines()
    for line in card:
        parsed_data.append(tuple(line.strip().split()))

    parsed_data.pop(0)  # removing header

    assert len(hubbard_data) == len(parsed_data)

    for array, parsed_array in zip(hubbard_data, parsed_data):
        assert array == parsed_array


@pytest.mark.parametrize(('parameters', 'values'), (
    (
        [[3, '1s', 3, '1s', 0.0, [0, 0, 0], 'U']],
        [
            ['Li', 'Co', 'O', 'O'],
            [[0, '1s', 0, '1s', 0.0, [0, 0, 0], 'U']],
        ],
    ),
    (
        [[1, '1s', 1, '1s', 0.0, [0, 0, 0], 'U'], [2, '2p', 2, '2p', 0.0, [0, 0, 0], 'U']],
        [
            ['O', 'O', 'Co', 'Li'],
            [[0, '1s', 0, '1s', 0.0, [0, 0, 0], 'U'], [1, '2p', 1, '2p', 0.0, [0, 0, 0], 'U']],
        ],
    ),
    (
        [[0, '3d', 0, '3d', 5.0, [0, 0, 0], 'U'], [3, '1s', 3, '1s', 0.0, [0, 0, 0], 'U']],
        [
            ['Li', 'Co', 'O', 'O'],
            [[0, '1s', 0, '1s', 0.0, [0, 0, 0], 'U'], [1, '3d', 1, '3d', 5.0, [0, 0, 0], 'U']],
        ],
    ),
    (
        [[1, '1s', 3, '2p', 0.0, [0, 0, 0], 'U']],
        [['O', 'O', 'Li', 'Co'], [[0, '1s', 2, '2p', 0.0, [0, 0, 0], 'U']]],
    ),
))
@pytest.mark.usefixtures('aiida_profile')
def test_reorder_atoms(generate_hubbard_structure, parameters, values):
    """Test the `reorder_atoms` method.

    .. note:: the test is strictly related to the sort logic implemented. If
        this was about to change, this test will most probably fail, and rearrangement
        of the input ``values`` should be performed.
    """
    # Sanity check - we must not override with default value
    projectors = 'atomic'
    formulation = 'liechtenstein'
    args = (parameters, projectors, formulation)
    hubbard_structure = generate_hubbard_structure(*args)
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)
    hubbard_utils.reorder_atoms()

    sites = hubbard_utils.hubbard_structure.sites
    for name, site in zip(values[0], sites):
        assert site.kind_name == name

    expected = Hubbard.from_list(values[1]).to_list()
    assert hubbard_utils.hubbard_structure.hubbard.projectors == projectors
    assert hubbard_utils.hubbard_structure.hubbard.formulation == formulation
    for param in hubbard_utils.hubbard_structure.hubbard.to_list():
        assert param in expected


@pytest.mark.usefixtures('aiida_profile')
def test_is_to_reorder(generate_hubbard_structure):
    """Test the `is_to_reorder` method."""
    parameters = [[1, '1s', 3, '2p', 0.0, [0, 0, 0], 'U']]
    hubbard_structure = generate_hubbard_structure(parameters)
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)
    assert hubbard_utils.is_to_reorder()

    parameters = [[0, '1s', 0, '2p', 0.0, [0, 0, 0], 'U']]
    hubbard_structure = generate_hubbard_structure(parameters)
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)
    assert not hubbard_utils.is_to_reorder()


@pytest.mark.usefixtures('aiida_profile')
def test_reorder_supercell_atoms(generate_hubbard_structure):
    """Test the `reorder_atoms` method with a supercell."""
    parameters = [
        (0, '3d', 0, '3d', 5.0, (0, 0, 0), 'U'),
        (0, '3d', 1, '2p', 5.0, (0, 0, 0), 'U'),
    ]
    hubbard_structure = generate_hubbard_structure(parameters=parameters)
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)

    pymatgen = hubbard_structure.get_pymatgen_structure()
    pymatgen.make_supercell([2, 2, 2])
    supercell = StructureData(pymatgen=pymatgen)
    hubbard_supercell = hubbard_utils.get_hubbard_for_supercell(supercell=supercell, thr=1e-5)
    supercell_utils = HubbardUtils(hubbard_supercell)
    supercell_utils.reorder_atoms()
    hubbard_supercell = supercell_utils.hubbard_structure

    for parameters in hubbard_supercell.hubbard.to_list():
        if parameters[0] == parameters[2]:
            assert hubbard_supercell.sites[parameters[0]].kind_name == 'Co'


@pytest.mark.usefixtures('aiida_profile')
def test_hubbard_for_supercell(generate_hubbard_structure):
    """Test the `get_hubbard_for_supercell` method.

    .. warning:: we only test for the ``U`` case, assuming that if U is transfered
        to the supercell correctly, any parameter will.
    """
    parameters = [[3, '1s', 3, '2s', 5.0, [0, 0, 0], 'U']]  # Li atom
    projectors = 'atomic'
    formulation = 'liechtenstein'
    args = (parameters, projectors, formulation)
    hubbard_structure = generate_hubbard_structure(*args)
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)

    pymatgen = hubbard_structure.get_pymatgen_structure()
    pymatgen.make_supercell([2, 2, 2])
    supercell = StructureData(pymatgen=pymatgen)

    hubbard_supercell = hubbard_utils.get_hubbard_for_supercell(supercell=supercell, thr=1e-5)

    assert hubbard_supercell.hubbard.projectors == projectors
    assert hubbard_supercell.hubbard.formulation == formulation
    for parameters in hubbard_supercell.hubbard.to_list():
        assert parameters[0] == parameters[2]
        assert parameters[1] == '1s'
        assert parameters[3] == '2s'
        assert hubbard_supercell.sites[parameters[0]].kind_name == 'Li'


@pytest.fixture
def get_non_trivial_hubbard_structure(filepath_tests):
    """Return a multi-coordination number `HubbardStructureData`."""

    def _get_non_trivial_hubbard_structure(name=None, use_kinds=False):
        """Return a multi-coordination number `HubbardStructureData`."""
        path = os.path.join(filepath_tests, 'fixtures', 'structures')

        if name is None:
            atoms = read(os.path.join(path, 'Fe3O4.cif'))

            if use_kinds:
                tags = [0] * atoms.get_global_number_of_atoms()
                for i in range(10):
                    tags[i] = 1
                tags[0] = 2
                for i in [26, 27, 28, 30, 31, 32]:
                    tags[i] = 1
                atoms.set_tags(tags)

            hubbard_structure = HubbardStructureData.from_structure(StructureData(ase=atoms))
            pymat = hubbard_structure.get_pymatgen_structure()

            tag_names = ['Fe']
            if use_kinds:
                tag_names = ['Fe', 'Fe1', 'Fe2']

            for i, site in enumerate(hubbard_structure.sites):
                if site.kind_name in tag_names:
                    hubbard_structure.append_hubbard_parameter(i, '3d', i, '3d', 5.0)
                    pymat_sites = pymat.get_sites_in_sphere(site.position, r=2.96)

                    for pymat_site in pymat_sites:
                        if pymat_site.position_atol != site.position:
                            # if pymat_site.specie.name != 'Fe':
                            translation = np.array(pymat_site.image, dtype=np.int64).tolist()
                            args = (i, '3d', int(pymat_site.index), '2p', 1.0, translation)
                            hubbard_structure.append_hubbard_parameter(*args)

        return hubbard_structure

    return _get_non_trivial_hubbard_structure


def test_get_interacting_pairs(get_non_trivial_hubbard_structure):
    """Test the `HubbardUtils.get_interacting_pairs` method."""
    hubbard_structure = get_non_trivial_hubbard_structure()
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)

    pairs = hubbard_utils.get_interacting_pairs()

    assert 'Fe' in pairs
    assert pairs['Fe'] == ['O']

    # using actual kinds
    hubbard_structure = get_non_trivial_hubbard_structure(use_kinds=True)
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)

    pairs = hubbard_utils.get_interacting_pairs()
    # print(hubbard_utils.get_hubbard_card())

    assert pairs == {'Fe': ['O', 'O1'], 'Fe1': ['O', 'O1'], 'Fe2': ['O', 'O1']}


def test_get_pairs_radius(get_non_trivial_hubbard_structure):
    """Test the `HubbardUtils.get_pairs_radius` method."""
    hubbard_structure = get_non_trivial_hubbard_structure()
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)

    rmin, rmax = hubbard_utils.get_pairs_radius(0, ['O'], 6)
    assert rmin < 2.2 < rmax

    rmin, rmax = hubbard_utils.get_pairs_radius(23, ['O'], 4)
    assert rmin < 2.2 < rmax

    # using actual kinds
    hubbard_structure = get_non_trivial_hubbard_structure(use_kinds=True)
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)

    rmin, rmax = hubbard_utils.get_pairs_radius(0, ['O', 'O1'], 6)
    assert rmin < 2.2 < rmax


def test_get_intersites_radius(get_non_trivial_hubbard_structure):
    """Test the `HubbardUtils.get_intersites_radius` method."""
    hubbard_structure = get_non_trivial_hubbard_structure()
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)

    radius = hubbard_utils.get_intersites_radius()

    assert abs(radius - 2.106) < 1.0e-5

    pymat = hubbard_structure.get_pymatgen_structure()
    for i in range(24):
        assert len(pymat.get_neighbors(pymat[i], r=radius)) in [4, 6]


def test_get_intersites_radius_five_nn(filepath_tests):
    """Test the `HubbardUtils.get_intersites_radius` method against LiMnTe (5 neighbours)."""
    path = os.path.join(filepath_tests, 'fixtures', 'structures', 'LMT.cif')
    atoms = read(path)

    hubbard_structure = HubbardStructureData.from_structure(StructureData(ase=atoms))
    hubbard_structure.initialize_onsites_hubbard('Mn', '3d')
    hubbard_structure.initialize_intersites_hubbard('Mn', '3d', 'Te', '2p')
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)

    radius = hubbard_utils.get_intersites_radius()

    pymat = hubbard_structure.get_pymatgen_structure()
    assert len(pymat.get_neighbors(pymat[0], r=radius)) == 5


def test_initialize_hubbard_parameters(get_non_trivial_hubbard_structure):
    """Test the `HubbardUtils.initialize_hubbard_parameters` method."""
    structure = get_non_trivial_hubbard_structure()
    hubbard_structure = initialize_hubbard_parameters(structure=structure, pairs={'Fe': ['3d', 5.0, 1e-8, {'O': '2p'}]})

    assert len(hubbard_structure.hubbard.parameters) == 8 * (4 + 1) + 16 * (6 + 1)


def test_initialize_hubbard_parameters_five_nn(filepath_tests):
    """Test the `HubbardUtils.initialize_hubbard_parameters` method against LiMnTe (5 neighbours)."""
    path = os.path.join(filepath_tests, 'fixtures', 'structures', 'LMT.cif')
    atoms = read(path)

    structure = StructureData(ase=atoms)

    hubbard_structure = initialize_hubbard_parameters(
        structure=structure, pairs={'Mn': ['3d', 5.0, 1e-8, {
            'Te': '4p'
        }]}
    )

    assert len(hubbard_structure.hubbard.parameters) == 6


@pytest.mark.parametrize(
    ('pairs', 'number_of_parameters'),
    (
        (
            {
                'Mn': ['3d', 5.0, 1e-8, {
                    'S': '3p'
                }],
                'Co': ['3d', 5.0, 1e-8, {
                    'S': '3p'
                }],
            },
            (1 + 6) + (1 + 4)  # 2 onsites, 6 + 4 intersites
        ),
        (
            {
                'Mn': ['3d', 5.0, 1e-8, {
                    'S': '3p',
                    'Co': '3d'
                }],
                'Co': ['3d', 5.0, 1e-8, {
                    'S': '3p',
                    'Mn': '3d'
                }],
            },
            (1 + 6 + 4) + (1 + 4 + 4)  # 2 onsites, 6 + 4 TM-S, 4 + 4 TM-TM
        ),
    )
)
def test_get_intersites_list(filepath_tests, pairs, number_of_parameters):
    """Test the `HubbardUtils.get_intersites_list` method against MnCoS."""
    path = os.path.join(filepath_tests, 'fixtures', 'structures', 'MnCoS.cif')
    atoms = read(path)

    structure = StructureData(ase=atoms)
    hubbard_structure = initialize_hubbard_parameters(structure=structure, pairs=pairs)

    assert len(hubbard_structure.hubbard.parameters) == number_of_parameters

    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)
    intersites = hubbard_utils.get_intersites_list()

    init_intersites = np.array(hubbard_structure.hubbard.to_list(), dtype='object')[:, [0, 2, 5]].tolist()

    assert len(intersites) == len(init_intersites)

    for intersite in intersites:
        assert intersite in init_intersites


def test_get_max_number_of_neighbours(filepath_tests):
    """Test the `HubbardUtiles.get_max_number_of_neighbours` method."""
    path = os.path.join(filepath_tests, 'fixtures', 'structures', 'MnCoS.cif')
    atoms = read(path)

    pairs = {
        'Mn': ['3d', 5.0, 1e-8, {
            'S': '3p',
        }],
        'Co': ['3d', 5.0, 1e-8, {
            'S': '3p',
        }],
    }
    structure = StructureData(ase=atoms)
    hubbard_structure = initialize_hubbard_parameters(structure=structure, pairs=pairs)
    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)

    assert hubbard_utils.get_max_number_of_neighbours() == 10


def test_max_number_of_neighbours(filepath_tests):
    """Test the `max_number_of_neighbours` method."""
    from aiida_quantumespresso.utils.hubbard import max_number_of_neighbours

    array = [[0, 1], [0, 1], [0, 0], [0, 1], [0, 2], [0, 2]]

    assert max_number_of_neighbours(array) == 5

    path = os.path.join(filepath_tests, 'fixtures', 'structures', 'MnCoS.cif')
    atoms = read(path)

    pairs = {
        'Mn': ['3d', 5.0, 1e-8, {
            'S': '3p',
            'Co': '3d'
        }],
        'Co': ['3d', 5.0, 1e-8, {
            'S': '3p',
            'Mn': '3d'
        }],
    }
    structure = StructureData(ase=atoms)
    hubbard_structure = initialize_hubbard_parameters(structure=structure, pairs=pairs)

    hubbard_utils = HubbardUtils(hubbard_structure=hubbard_structure)
    intersites = np.array(hubbard_utils.get_intersites_list(), dtype='object')[:, [0, 1]]

    assert max_number_of_neighbours(intersites) == 10
