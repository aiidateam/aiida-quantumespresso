# -*- coding: utf-8 -*-
"""Tests for the `create_kpoints_from_distance` calculation function."""

from aiida import orm

from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import create_kpoints_from_distance


def test_3d_cubic(generate_structure):
    """Test `create_kpoints_from_distance` calculation function.

    Case: cubic cell with all sides equal and all periodic.
    Expected result: isotropic kpoint mesh.
    """
    structure = generate_structure('silicon-kinds')

    kpoints = create_kpoints_from_distance(structure, orm.Float(0.3), orm.Bool(False)).get_kpoints_mesh()[0]

    assert kpoints[0] == kpoints[1] == kpoints[2]


def test_2d_vacuum(generate_structure):
    """Test `create_kpoints_from_distance` calculation function.

    Case: periodicity only in two directions, third direction has vacuum (long lattice vector).
    Expected result: anisotropic kpoint mesh with third direction set to 1.
    """
    structure = generate_structure('2D-xy-arsenic')

    kpoints = create_kpoints_from_distance(structure, orm.Float(0.3), orm.Bool(False)).get_kpoints_mesh()[0]

    assert kpoints[0] == kpoints[1]
    assert kpoints[2] == 1


def test_2d_cubic(generate_structure):
    """Test `create_kpoints_from_distance` calculation function.

    Case: cubic cell with all sides equal and periodicity is 2D.
    Expected result: non periodic direction is set to 1.
    """
    structure = generate_structure('silicon-kinds')
    structure.pbc = (True, True, False)

    kpoints = create_kpoints_from_distance(structure, orm.Float(0.3), orm.Bool(False)).get_kpoints_mesh()[0]

    assert kpoints[0] == kpoints[1]
    assert kpoints[2] == 1
