# -*- coding: utf-8 -*-
"""Create a new magnetic configuration from the given structure based on the desired magnetic moments."""
from aiida.engine import calcfunction
from aiida.orm import Float
import numpy


@calcfunction
def create_magnetic_configuration(
    structure, magnetic_moment_per_site, atol=lambda: Float(0.5), ztol=lambda: Float(0.05)
):
    """Create a new magnetic configuration from the given structure based on a list of magnetic moments per site.

    To create the new list of kinds, the algorithm loops over all the elements in the structure and makes a list of the
    sites with that element and their corresponding magnetic moment. Next, it splits this list in three lists:

    * Zero magnetic moments: Any site that has an absolute magnetic moment lower than ``ztol``
    * Positive magnetic moments
    * Negative magnetic moments

    The algorithm then sorts the positive and negative lists from large to small absolute value, and loops over each of
    list. New magnetic kinds will be created when the absolute difference between the magnetic moment of the current
    kind and the site exceeds ``atol``.

    The positive and negative magnetic moments are handled separately to avoid assigning two sites with opposite signs
    in their magnetic moment to the same kind and make sure that each kind has the correct magnetic moment, i.e. the
    largest magnetic moment in absolute value of the sites corresponding to that kind.

    .. important:: the function currently does not support alloys.

    :param structure: a `StructureData` instance.
    :param magnetic_moment_per_site: list of magnetic moments for each site in the structure.
    :param atol: the absolute tolerance on determining if two sites have the same magnetic moment.
    :param ztol: threshold for considering a kind to have non-zero magnetic moment.
    """
    # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    import string

    from aiida.orm import Dict, StructureData

    if structure.is_alloy:
        raise ValueError('Alloys are currently not supported.')

    atol = atol.value
    rtol = 0  # Relative tolerance used in the ``numpy.is_close()`` calls.
    ztol = ztol.value

    new_structure = StructureData(cell=structure.cell, pbc=structure.pbc)
    magnetic_configuration = {}

    for element in structure.get_symbols_set():

        # Filter the sites and magnetic moments on the site element
        element_sites, element_magnetic_moments = zip(
            *[(site, magnetic_moment)
              for site, magnetic_moment in zip(structure.sites, magnetic_moment_per_site)
              if site.kind_name.rstrip(string.digits) == element]
        )

        # Split the sites and their magnetic moments by sign to filter out the sites with magnetic moment lower than
        # `ztol`and deal with the positive and negative magnetic moment separately. This is important to avoid assigning
        # two sites with opposite signs to the same kind and make sure that each kind has the correct magnetic moment,
        # i.e. the largest magnetic moment in absolute value of the sites corresponding to that kind.
        zero_sites = []
        pos_sites = []
        neg_sites = []

        for site, magnetic_moment in zip(element_sites, element_magnetic_moments):

            if abs(magnetic_moment) <= ztol:
                zero_sites.append((site, 0))
            elif magnetic_moment > 0:
                pos_sites.append((site, magnetic_moment))
            else:
                neg_sites.append((site, magnetic_moment))

        kind_index = -1
        kind_names = []
        kind_sites = []
        kind_magnetic_moments = {}

        for site_list in (zero_sites, pos_sites, neg_sites):

            if not site_list:
                continue

            # Sort the site list in order to build the kind lists from large to small absolute magnetic moment.
            site_list = sorted(site_list, key=lambda x: abs(x[1]), reverse=True)

            sites, magnetic_moments = zip(*site_list)

            kind_index += 1
            current_kind_name = f'{element}{kind_index}'
            kind_sites.append(sites[0])
            kind_names.append(current_kind_name)
            kind_magnetic_moments[current_kind_name] = magnetic_moments[0]

            for site, magnetic_moment in zip(sites[1:], magnetic_moments[1:]):

                if not numpy.isclose(magnetic_moment, kind_magnetic_moments[current_kind_name], rtol, atol):
                    kind_index += 1
                    current_kind_name = f'{element}{kind_index}'
                    kind_magnetic_moments[current_kind_name] = magnetic_moment

                kind_sites.append(site)
                kind_names.append(current_kind_name)

        # In case there is only a single kind for the element, remove the 0 kind index
        if current_kind_name == f'{element}0':
            kind_names = len(element_magnetic_moments) * [element]
            kind_magnetic_moments = {element: kind_magnetic_moments[current_kind_name]}

        magnetic_configuration.update(kind_magnetic_moments)

        for name, site in zip(kind_names, kind_sites):
            new_structure.append_atom(
                name=name,
                symbols=(element,),
                weights=(1.0,),
                position=site.position,
            )

    return {'structure': new_structure, 'magnetic_moments': Dict(dict=magnetic_configuration)}
