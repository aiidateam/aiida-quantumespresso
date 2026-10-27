"""Tests for the `rattle_structure` calculation function."""

import numpy as np
from aiida.orm import Float

from aiida_quantumespresso.calculations.functions.rattle_structure import rattle_structure
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData


def test_rattle_structure(generate_structure):
    """Test the `rattle_structure` calculation function with a plain ``StructureData``."""
    structure = generate_structure('silicon')
    rattled = rattle_structure(structure, Float(0.01))

    assert rattled.cell == structure.cell
    assert [site.kind_name for site in rattled.sites] == [site.kind_name for site in structure.sites]

    positions = np.array([site.position for site in structure.sites])
    rattled_positions = np.array([site.position for site in rattled.sites])
    assert not np.allclose(positions, rattled_positions)
    assert np.allclose(positions, rattled_positions, atol=0.1)


def test_rattle_structure_hubbard(generate_structure):
    """Test that rattling a ``HubbardStructureData`` preserves the node type and the Hubbard parameters."""
    hubbard_structure = HubbardStructureData.from_structure(generate_structure('silicon'))
    hubbard_structure.initialize_onsites_hubbard('Si', '3p', 5.0)

    rattled = rattle_structure(hubbard_structure, Float(0.01))

    assert isinstance(rattled, HubbardStructureData)
    assert rattled.hubbard == hubbard_structure.hubbard
