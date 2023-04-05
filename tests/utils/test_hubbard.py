# -*- coding: utf-8 -*-
"""Tests for the :mod:`utils.hubbard` module."""
# pylint: disable=redefined-outer-name
import os

import pytest

from aiida_quantumespresso.common.hubbard import Hubbard
from aiida_quantumespresso.utils.hubbard import HubbardUtils


@pytest.fixture
def generate_hubbard_structure():
    """Return a `HubbardStructureData` instance."""

    def _generate_hubbard_structure(parameters=None, projectors=None, formulation=None):
        from aiida.orm import StructureData

        from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

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

    def _generate_hubbard_utils():
        return HubbardUtils(hubbard_structure=generate_hubbard_structure())

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
    from aiida.orm import StructureData
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
    from aiida.orm import StructureData
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
