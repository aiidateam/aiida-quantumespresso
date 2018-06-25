# -*- coding: utf-8 -*-
from aiida.work.workfunctions import workfunction


@workfunction
def create_kpoints_from_distance(structure, distance, force_parity):
    """
    Generate a uniformly spaced kpoint mesh for a given structure where the spacing between kpoints in reciprocal
    space is guaranteed to be at least the defined distance.

    :param structure: the StructureData to which the mesh should apply
    :param distance: a Float with the desired distance between kpoints in reciprocal space
    :param force_parity: a Bool to specify whether the generated mesh should maintain parity
    :returns: a KpointsData with the generated mesh
    """
    from aiida.orm.data.array.kpoints import KpointsData

    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    kpoints.set_kpoints_mesh_from_density(distance.value, force_parity=force_parity.value)

    return kpoints
