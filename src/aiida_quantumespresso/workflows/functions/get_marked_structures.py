# -*- coding: utf-8 -*-
"""CalcFunction to create structures with a marked atom for each site in a list."""
import warnings

from aiida import orm
from aiida.common import ValidationError
from aiida.engine import calcfunction
from aiida.orm.nodes.data.structure import Kind, Site, StructureData

warnings.warn(
    'This module is deprecated and will be removed soon as part of migrating XAS and XPS workflows to a new repository.'
    '\nThe new repository can be found at: https://github.com/aiidaplugins/aiida-qe-xspec.', FutureWarning
)


@calcfunction
def get_marked_structures(structure, atoms_list, marker='X'):
    """Read a StructureData object and return structures for XPS calculations.

    :param atoms_list:  the atoms_list of atoms to be marked.
    :param marker: a Str node defining the name of the marked atom Kind. Default is 'X'.
    :returns: StructureData objects for the generated structure.
    """
    marker = marker.value
    elements_present = [kind.symbol for kind in structure.kinds]
    if marker in elements_present:
        raise ValidationError(
            f'The marker ("{marker}") should not match an existing Kind in '
            f'the input structure ({elements_present}.'
        )

    output_params = {}
    result = {}

    for index in atoms_list.get_list():
        marked_structure = StructureData()
        kinds = {kind.name: kind for kind in structure.kinds}
        marked_structure.set_cell(structure.cell)

        for i, site in enumerate(structure.sites):
            if i == index:
                marked_kind = Kind(name=marker, symbols=site.kind_name)
                marked_site = Site(kind_name=marked_kind.name, position=site.position)
                marked_structure.append_kind(marked_kind)
                marked_structure.append_site(marked_site)
                output_params[f'site_{index}'] = {'symbol': site.kind_name, 'multiplicity': 1}
            else:
                if site.kind_name not in [kind.name for kind in marked_structure.kinds]:
                    marked_structure.append_kind(kinds[site.kind_name])
                new_site = Site(kind_name=site.kind_name, position=site.position)
                marked_structure.append_site(new_site)
        result[f'site_{index}'] = marked_structure

    result['output_parameters'] = orm.Dict(dict=output_params)

    return result
