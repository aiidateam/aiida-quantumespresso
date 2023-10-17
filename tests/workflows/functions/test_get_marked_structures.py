# -*- coding: utf-8 -*-
"""Tests for the `get_marked_structure` class."""


def test_get_marked_structure():
    """Test the get_marked_structure function."""
    from aiida.orm import List, StructureData
    from ase.build import molecule

    from aiida_quantumespresso.workflows.functions.get_marked_structures import get_marked_structures

    mol = molecule('CH3CH2OH')
    mol.center(vacuum=2.0)
    structure = StructureData(ase=mol)
    indices = List(list=[0, 1, 2])
    output = get_marked_structures(structure, indices)
    assert len(output) == 4
    assert output['site_0'].get_site_kindnames() == ['X', 'C', 'O', 'H', 'H', 'H', 'H', 'H', 'H']
    assert output['site_1'].get_site_kindnames() == ['C', 'X', 'O', 'H', 'H', 'H', 'H', 'H', 'H']
