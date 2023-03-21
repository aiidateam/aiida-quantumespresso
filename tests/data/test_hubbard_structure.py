# -*- coding: utf-8 -*-
"""Tests for the :mod:`data.hubbard_structure` module."""
import pytest
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

@pytest.fixture
def generate_hubbard():
    """Return a `Hubbard` instance."""
    
    def _generate_hubbard():
        from aiida_quantumespresso.common.hubbard import Hubbard
        return Hubbard.from_list([[0,'1s',0,'1s',5.0,[0,0,0],'Dudarev-Ueff']])
    
    return _generate_hubbard

@pytest.fixture
def generate_hubbard_structure(generate_structure):
    """Return a `HubbardStructureData` instance."""
    
    def _generate_hubbard_structure():
        from aiida_quantumespresso.common.hubbard import Hubbard
        structure = generate_structure('silicon-kinds')
        hp_list = [[0,'1s',0,'1s',5.0,[0,0,0],'Dudarev-Ueff']]
        hubbard = Hubbard.from_list(hp_list)
        return HubbardStructureData.from_structure(structure=structure, hubbard=hubbard)
    
    return _generate_hubbard_structure

@pytest.mark.usefixtures('aiida_profile')
def test_valid_init(generate_hubbard):
    """Test the constructor."""
    cell = [[1,0,0],[0,1,0],[0,0,1]]
    sites = [['Si', 'Si0', [0,0,0]]]
    hubbard = generate_hubbard()
    hs = HubbardStructureData(cell=cell, sites=sites, hubbard=hubbard)
    assert hs.cell == cell
    assert hs.kinds[0].symbol == sites[0][0]
    assert hs.sites[0].kind_name == sites[0][1]
    assert list(hs.sites[0].position) == sites[0][2]
    assert hs.hubbard == hubbard

@pytest.mark.usefixtures('aiida_profile')
def test_from_structure(generate_structure, generate_hubbard):
    """Test the `from_structure` method."""
    structure = generate_structure()
    hubbard = generate_hubbard()
    hs = HubbardStructureData.from_structure(structure=structure, hubbard=hubbard)
    assert hs.cell == structure.cell
    assert len(hs.sites) == len(structure.sites)
    assert hs.hubbard == hubbard # sanity check 

    structure = generate_structure('silicon-kinds')
    hs = HubbardStructureData.from_structure(structure=structure)
    assert hs.get_site_kindnames() == structure.get_site_kindnames()
    assert len(hs.kinds) == 2

@pytest.mark.usefixtures('aiida_profile')
def test_append_hubbard_parameters(generate_hubbard_structure):
    """Test the `append_hubbard_parameters` method."""
    from aiida_quantumespresso.common.hubbard import HubbardParameters
    hs = generate_hubbard_structure()
    args = [0,'1s',1,'1s', 5.0, [0,0,0], 'Dudarev-U']
    hs.append_hubbard_parameter(*args)
    hp = HubbardParameters.from_list(args)
    assert len(hs.hubbard.hubbard_parameters) == 2
    assert hp == hs.hubbard.hubbard_parameters[1]

@pytest.mark.usefixtures('aiida_profile')
def test_pop_hubbard_parameter(generate_hubbard_structure):
    """Test the `pop_hubbard_parameter` method."""
    hs = generate_hubbard_structure()
    hs.pop_hubbard_parameter(0)
    assert len(hs.hubbard.hubbard_parameters) == 0

@pytest.mark.usefixtures('aiida_profile')
def test_clear_hubbard_parameters(generate_hubbard_structure):
    """Test the `clear_hubbard_parameters` method."""
    hs = generate_hubbard_structure()
    hs.clear_hubbard_parameters()
    assert len(hs.hubbard.hubbard_parameters) == 0

@pytest.mark.usefixtures('aiida_profile')
def test_is_storable(generate_hubbard_structure):
    """Test the storing does not throw errors."""
    hs = generate_hubbard_structure()
    hs.store()
    assert hs.is_stored

@pytest.mark.usefixtures('aiida_profile')
def test_initialize_intersites_hubbard(generate_hubbard_structure):
    """Test the `initialize_intersites_hubbard` method."""
    hs = generate_hubbard_structure()
    hs.initialize_intersites_hubbard('Si','1s','Si','2s',0,'Dudarev-Ueff',False)
    
    # !WARNING! This is not the expected behavior, as we would like it to initialize
    #           intersites among first neighbours. The method was designed for different
    #           interacting species. We may want to improve it for this special cases,
    #           although it is still rather futuristic.
    assert hs.hubbard.hubbard_parameters[1].to_list() == [0,'1s',0,'2s',0.0,[0,0,0],'Dudarev-Ueff'] 
    
    hs.clear_hubbard_parameters()
    hs.initialize_intersites_hubbard('Si0','1s','Si1','2s',0.0,'Dudarev-Ueff')
    assert [0,'1s',1,'2s',0,[-1,0,0],'Dudarev-Ueff']  in hs.hubbard.to_list()
    assert len(hs.hubbard.hubbard_parameters) == 1

    with pytest.raises(ValueError):
        hs.initialize_intersites_hubbard('Mg','1s','Si1','2s',0,'Dudarev-Ueff')

@pytest.mark.usefixtures('aiida_profile')
def test_initialize_onsites_hubbard(generate_hubbard_structure):
    """Test the `initialize_onsites_hubbard` method."""
    hs = generate_hubbard_structure()

    hs.clear_hubbard_parameters()
    hs.initialize_onsites_hubbard('Si','1s',0.0,'Dudarev-Ueff',False)
   
    assert [0,'1s',0,'1s',0,[0,0,0],'Dudarev-Ueff'] in hs.hubbard.to_list()
    assert len(hs.hubbard.hubbard_parameters) == 2

    hs.clear_hubbard_parameters()
    hs.initialize_onsites_hubbard('Si0','1s',0.0,'Dudarev-Ueff',True)

    assert len(hs.hubbard.hubbard_parameters) == 1

@pytest.mark.usefixtures('aiida_profile')
def test_get_one_kind_index(generate_hubbard_structure):
    """Test the `_get_one_kind_index` method."""
    hs = generate_hubbard_structure()
    hs._get_one_kind_index('Si0') == [0]
    hs._get_one_kind_index('Si1') == [1]

@pytest.mark.usefixtures('aiida_profile')
def test_get_symbol_indecis(generate_hubbard_structure):
    """Test the `_get_symbol_indecis` method."""
    hs = generate_hubbard_structure()
    hs._get_symbol_indecis('Si') == [0]

