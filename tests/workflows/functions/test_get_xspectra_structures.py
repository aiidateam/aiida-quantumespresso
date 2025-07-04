# -*- coding: utf-8 -*-
"""Tests for the `get_marked_structure` class."""
from aiida.orm import Bool, Dict, Int
import pytest

from aiida_quantumespresso.utils.hubbard import HubbardStructureData, HubbardUtils
from aiida_quantumespresso.workflows.functions.get_xspectra_structures import get_xspectra_structures


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


def test_base(generate_structure):
    """Test the basic operation of get_xspectra_structures."""

    c_si = generate_structure('silicon')
    spglib_settings = Dict({'symprec': 1.0e-3})
    inputs = {'structure': c_si, 'spglib_settings': spglib_settings}
    result = get_xspectra_structures(**inputs)
    assert len(result) == 4
    out_params = result['output_parameters'].get_dict()
    assert out_params['spacegroup_number'] == 227
    assert out_params['supercell_num_sites'] == 64
    assert len(out_params['equivalent_sites_data']) == 1


def test_use_element_types(generate_structure):
    """Test the CF's `use_element_types` flag."""

    c_si = generate_structure('silicon')
    c_si_kinds = generate_structure('silicon-kinds')
    spglib_settings = Dict({'symprec': 1.0e-3})
    inputs_bare = {'structure': c_si, 'spglib_settings': spglib_settings}
    inputs_kinds = {
        'structure': c_si_kinds,
        'use_element_types': Bool(False),
        'spglib_options': spglib_settings,
        'standardize_structure': Bool(False)
    }

    result_bare = get_xspectra_structures(**inputs_bare)
    result_kinds = get_xspectra_structures(**inputs_kinds)

    inputs_element_types = {
        'structure': c_si_kinds,
        'spglib_options': spglib_settings,
        'standardize_structure': Bool(False),
    }
    result_element_types = get_xspectra_structures(**inputs_element_types)

    assert 'site_1_Si' in result_kinds
    assert 'site_1_Si' not in result_element_types
    assert 'site_1_Si' not in result_bare


def test_hubbard(generate_structure):
    """Test that the CalcFunction will pass Hubbard parameters to the output structures.

    The intent here is to confirm that simply using the `initialize_` methods to get
    hubbard parameters will propogate to the resulting supercells and (crucially)
    generate the correct hubbard card.
    """

    c_si = generate_structure()
    c_si_hub = HubbardStructureData.from_structure(c_si)
    c_si_hub.initialize_onsites_hubbard('Si', '1s', 0.0, 'Ueff', False)

    inputs = {'structure': c_si_hub, 'standardize_structure': Bool(False)}
    result = get_xspectra_structures(**inputs)

    marked = result['site_0_Si']
    utils_marked = HubbardUtils(marked)
    hub_card_lines = [i.strip() for i in utils_marked.get_hubbard_card().splitlines()]
    out_params = result['output_parameters'].get_dict()

    assert out_params['supercell_num_sites'] == 54
    assert 'U\tX-1s\t0.0' in hub_card_lines
    assert 'U\tSi-1s\t0.0' in hub_card_lines


def test_symmetry_input(generate_structure):
    """Test the basic operation of get_xspectra_structures."""

    c_si = generate_structure('silicon')
    sites_data = {'site_0': {'symbol': 'Si', 'site_index': 1, 'multiplicity': 8}}
    inputs = {
        'structure': c_si,
        'equivalent_sites_data': Dict(sites_data),
        'spacegroup_number': Int(220),
        'parse_symmetry': Bool(False),
    }

    # Test also to confirm at leaving out `parse_symmetry` should cause
    # the CF to ignore the provided symmetry data.
    inputs_no_parse_symmetry = {
        'structure': c_si,
        'equivalent_sites_data': Dict(sites_data),
        'spacegroup_number': Int(220),
    }

    result = get_xspectra_structures(**inputs)
    result_no_parse_symmetry = get_xspectra_structures(**inputs_no_parse_symmetry)

    out_params = result['output_parameters'].get_dict()
    out_params_no_parse_symmetry = result_no_parse_symmetry['output_parameters'].get_dict()

    assert out_params['spacegroup_number'] == 220
    assert out_params['symmetry_parsed'] is False
    assert out_params['equivalent_sites_data'] == sites_data
    assert out_params_no_parse_symmetry['spacegroup_number'] == 227
    assert out_params_no_parse_symmetry['symmetry_parsed'] is True
    assert out_params_no_parse_symmetry['equivalent_sites_data'] != sites_data
