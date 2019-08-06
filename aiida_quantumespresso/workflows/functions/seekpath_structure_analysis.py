# -*- coding: utf-8 -*-
"""Calcfunction to primitivize a structure and return high symmetry k-point path through its Brillouin zone."""
from __future__ import absolute_import
from aiida.engine import calcfunction


@calcfunction
def seekpath_structure_analysis(structure, parameters):
    """Primitivize the structure with SeeKpath and return high symmetry k-point path through its Brillouin zone.

    This calcfunction will take a structure and pass it through SeeKpath to get the normalized primitive cell and the
    path of high symmetry k-points through its Brillouin zone. Note that the returned primitive cell may differ from the
    original structure in which case the k-points are only congruent with the primitive cell.
    """
    from aiida.tools import get_explicit_kpoints_path
    return get_explicit_kpoints_path(structure, **parameters.get_dict())
