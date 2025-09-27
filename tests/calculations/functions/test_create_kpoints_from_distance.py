# -*- coding: utf-8 -*-
"""Tests for the `create_kpoints_from_distance` calculation function."""
from aiida.orm import Bool, Float

from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import create_kpoints_from_distance


def test_kpoints_00(generate_structure):
    """Test `create_kpoints_from_distance` calculation function.

    Case: cubic cell with all sides equal and all periodic.
    Expected result: isotropic kpoint mesh.
    """
    structure = generate_structure('silicon-kinds')

    distance = Float(0.3)
    force_parity = Bool(False)

    kpoints = create_kpoints_from_distance(structure, distance, force_parity).get_kpoints_mesh()[0]

    assert kpoints[0] == kpoints[1] == kpoints[2]


def test_kpoints_01(generate_structure):
    """Test `create_kpoints_from_distance` calculation function.

    Case: periodicity only in two directions.
    Expected result: anisotropic kpoint mesh with one of the directions set to 1.
    """
    structure = generate_structure('2D-xy-arsenic')

    distance = Float(0.3)
    force_parity = Bool(False)

    kpoints = create_kpoints_from_distance(structure, distance, force_parity).get_kpoints_mesh()[0]

    assert kpoints[0] == kpoints[1]
    assert kpoints[2] == 1


def test_kpoints_02(generate_structure):
    """Test `create_kpoints_from_distance` calculation function.

    Case: cubic cell with all sides equal and periodicity is 2D.
    Expected result: non periodic direction is set to 1.
    """
    structure = generate_structure('silicon-kinds')
    structure.pbc = (True, True, False)

    distance = Float(0.3)
    force_parity = Bool(False)

    kpoints = create_kpoints_from_distance(structure, distance, force_parity).get_kpoints_mesh()[0]

    assert kpoints[0] == kpoints[1]
    assert kpoints[2] == 1
