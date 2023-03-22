# -*- coding: utf-8 -*-
"""Tests for the :mod:`utils.hubbard` module."""
import os
import pytest


@pytest.fixture
def generate_hubbard_structure():
    """Return a `HubbardStructureData` instance."""
    
    def _generate_hubbard_structure():
        from aiida.orm import StructureData
        from aiida_quantumespresso.common.hubbard import Hubbard
        from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

        a, b, c, d = 1.40803, 0.81293, 4.68453, 1.62585
        cell = [[a, -b, c], [0.0, d, c], [-a, -b, c]]
        positions = [[0, 0, 0], [0, 0, 3.6608], [0, 0, 10.392], [0, 0, 7.0268]]
        symbols = ['Co', 'O', 'O', 'Li']
        structure = StructureData(cell=cell)
        for position, symbol in zip(positions, symbols):
            structure.append_atom(position=position, symbols=symbol)

        hp_list = [[0,'3d',0,'3d',5.0,[0,0,0],'U']]
        hp_list = [[0,'3d',1,'2p',1.0,[0,0,0],'V']]
        hubbard = Hubbard.from_list(hp_list)

        return HubbardStructureData.from_structure(structure=structure, hubbard=hubbard)
    
    return _generate_hubbard_structure


@pytest.fixture
def generate_hubbard_utils(generate_hubbard_structure):
    """Return a `HubbardUtils` instance."""
    
    def _generate_hubbard_utils():
        from aiida_quantumespresso.utils.hubbard import HubbardUtils

        return HubbardUtils(hubbard_structure=generate_hubbard_structure())
    
    return _generate_hubbard_utils

@pytest.fixture
def filepath_hubbard(filepath_tests):
    """Return the absolute filepath to the directory containing the file `fixtures`."""
    return os.path.join(filepath_tests, 'utils/fixtures/hubbard')


@pytest.mark.usefixtures('aiida_profile')
def test_valid_init(generate_hubbard_utils, generate_hubbard_structure):
    """Test the constructor."""
    hubbard_utils = generate_hubbard_utils()
    assert hubbard_utils.hubbard_structure.hubbard == generate_hubbard_structure().hubbard


@pytest.mark.parametrize(
    ('filename','projectors'),
    (
        ("HUBBARD.dat", "ortho-atomic"),
        ("HUBBARD_2.dat", "atomic"),
    )
)
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
                hubbard_data.append(list(line.strip().split()))

    hubbard_data.pop(0) # removing header
    
    parsed_data = []
    card = hubbard_utils.get_hubbard_card()
    card = card.splitlines()
    for line in card:
        parsed_data.append(list(line.strip().split()))
    
    parsed_data.pop(0) # removing header
    
    assert len(hubbard_data)==len(parsed_data)
    
    for array, parsed_array in zip(hubbard_data, parsed_data):
        assert array==parsed_array



