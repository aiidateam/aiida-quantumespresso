# -*- coding: utf-8 -*-
"""Tests for the `seekpath_structure_analysis` function for HubbbardStructureData."""
import pytest

from aiida_quantumespresso.calculations.functions.seekpath_structure_analysis import seekpath_structure_analysis
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData


# pylint: disable=W0621
def test_seekpath_analysis(data_regression):
    """Test the `seekpath_structure_analysis` calculation function for HubbardStructureData."""
    cell = [[5.43, 0.0, 0.0], [0.0, 5.43, 0.0], [0.0, 0.0, 5.43]]
    sites = (
        ('Co', 'Co0', (0.0, 0.0, 0.0)),
        ('Co', 'Co0', (0.0, 2.715, 2.715)),
        ('Co', 'Co1', (2.715, 0.0, 2.715)),
        ('Co', 'Co1', (2.715, 2.715, 0.0)),
    )
    orig_structure = HubbardStructureData(cell=cell, sites=sites)

    for hubbard_parameter in [
        ('Co0', '3d', 4.0, 'Ueff', True),
        ('Co1', '3d', 3.0, 'U', True),
    ]:
        orig_structure.initialize_onsites_hubbard(*hubbard_parameter)

    result = seekpath_structure_analysis(orig_structure)

    prim_structure = result['primitive_structure']
    conv_structure = result['conv_structure']

    assert isinstance(prim_structure, HubbardStructureData), 'Primitive structure should be a HubbardStructureData'
    assert isinstance(conv_structure, HubbardStructureData), 'Conventional structure should be a HubbardStructureData'

    assert prim_structure.hubbard.parameters != orig_structure.hubbard.parameters, \
        'Primitive parameters should be different'
    assert len(prim_structure.hubbard.parameters) == len(orig_structure.hubbard.parameters), \
        'Primitive parameters should have the same length as original parameters'
    assert all(
        prim_param.atom_manifold == orig_param.atom_manifold
        for prim_param, orig_param in zip(prim_structure.hubbard.parameters, orig_structure.hubbard.parameters)
    ), 'Primitive cell parameter atom manifolds should match the original'

    data_regression.check({
        'primitive': {
            'cell': prim_structure.cell,
            'kinds': prim_structure.get_site_kindnames(),
            'positions': [site.position for site in prim_structure.sites],
            'hubbard': prim_structure.hubbard.to_list(),
        },
        'conventional': {
            'cell': conv_structure.cell,
            'kinds': conv_structure.get_site_kindnames(),
            'positions': [site.position for site in conv_structure.sites],
            'hubbard': conv_structure.hubbard.to_list(),
        }
    })


# pylint: disable=W0621
def test_seekpath_analysis_intersite(generate_structure):
    """Test that the `seekpath_structure_analysis` with intersite hubbard corrections fails."""
    orig_structure = HubbardStructureData.from_structure(generate_structure('silicon-kinds'))
    orig_structure.initialize_intersites_hubbard('Si0', '2p', 'Si1', '2p', 4.0, 'V', True)

    with pytest.raises(NotImplementedError, match='Intersite Hubbard parameters are not yet supported.'):
        seekpath_structure_analysis(orig_structure)
