# -*- coding: utf-8 -*-
from __future__ import absolute_import
from aiida.engine import calcfunction


@calcfunction
def seekpath_structure_analysis(structure, parameters):
    """
    This calcfunction will take a structure and pass it through SeeKpath to get the
    primitive cell and the path of high symmetry k-points through its Brillouin zone.
    Note that the returned primitive cell may differ from the original structure in
    which case the k-points are only congruent with the primitive cell.
    """
    from aiida.tools import get_explicit_kpoints_path
    return get_explicit_kpoints_path(structure, **parameters.get_dict())
