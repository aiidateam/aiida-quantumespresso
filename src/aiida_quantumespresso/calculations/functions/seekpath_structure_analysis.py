# -*- coding: utf-8 -*-
"""Calcfunction to primitivize a structure and return high symmetry k-point path through its Brillouin zone."""
from aiida.engine import calcfunction
from aiida.orm import Data

from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData


@calcfunction
def seekpath_structure_analysis(structure, **kwargs):
    """Primitivize the structure with SeeKpath and generate the high symmetry k-point path through its Brillouin zone.

    This calcfunction will take a structure and pass it through SeeKpath to get the normalized primitive cell and the
    path of high symmetry k-points through its Brillouin zone. Note that the returned primitive cell may differ from the
    original structure in which case the k-points are only congruent with the primitive cell.

    The keyword arguments can be used to specify various Seekpath parameters, such as:

        with_time_reversal: True
        reference_distance: 0.025
        recipe: 'hpkot'
        threshold: 1e-07
        symprec: 1e-05
        angle_tolerance: -1.0

    Note that exact parameters that are available and their defaults will depend on your Seekpath version.
    """
    from aiida.tools import get_explicit_kpoints_path

    # All keyword arugments should be `Data` node instances of base type and so should have the `.value` attribute
    unwrapped_kwargs = {key: node.value for key, node in kwargs.items() if isinstance(node, Data)}

    result = get_explicit_kpoints_path(structure, **unwrapped_kwargs)

    if isinstance(structure, HubbardStructureData):
        result['primitive_structure'] = update_structure_with_hubbard(result['primitive_structure'], structure)
        result['conv_structure'] = update_structure_with_hubbard(result['conv_structure'], structure)

    return result


def update_structure_with_hubbard(structure, orig_structure):
    """Update the structure based on Hubbard parameters if the input structure is a HubbardStructureData."""
    from aiida_quantumespresso.utils.hubbard import is_intersite_hubbard

    hubbard_structure = HubbardStructureData.from_structure(structure)

    if is_intersite_hubbard(orig_structure.hubbard):
        raise NotImplementedError('Intersite Hubbard parameters are not yet supported.')

    for parameter in orig_structure.hubbard.parameters:
        hubbard_structure.initialize_onsites_hubbard(
            atom_name=orig_structure.sites[parameter.atom_index].kind_name,
            atom_manifold=parameter.atom_manifold,
            value=parameter.value,
            hubbard_type=parameter.hubbard_type,
            use_kinds=True,
        )

    return hubbard_structure
