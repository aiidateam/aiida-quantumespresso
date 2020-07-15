# -*- coding: utf-8 -*-
"""Module with utitlies for the CLI to generate default values."""


def get_structure():
    """Return a `StructureData` representing bulk silicon.

    The database will first be queried for the existence of a bulk silicon crystal. If this is not the case, one is
    created and stored. This function should be used as a default for CLI options that require a `StructureData` node.
    This way new users can launch the command without having to construct or import a structure first. This is the
    reason that we hardcode a bulk silicon crystal to be returned. More flexibility is not required for this purpose.

    :return: a `StructureData` representing bulk silicon
    """
    from ase.spacegroup import crystal
    from aiida.orm import QueryBuilder, StructureData

    # Filters that will match any elemental Silicon structure with 2 or less sites in total
    filters = {
        'attributes.sites': {
            'of_length': 2
        },
        'attributes.kinds': {
            'of_length': 1
        },
        'attributes.kinds.0.symbols.0': 'Si'
    }

    builder = QueryBuilder().append(StructureData, filters=filters)
    results = builder.first()

    if not results:
        alat = 5.43
        ase_structure = crystal(
            'Si',
            [(0, 0, 0)],
            spacegroup=227,
            cellpar=[alat, alat, alat, 90, 90, 90],
            primitive_cell=True,
        )
        structure = StructureData(ase=ase_structure)
        structure.store()
    else:
        structure = results[0]

    return structure.uuid
