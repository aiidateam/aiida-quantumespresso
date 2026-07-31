"""Calculation function to apply random displacements to the atomic positions of a structure."""

from aiida import orm
from aiida.engine import calcfunction


@calcfunction
def rattle_structure(structure, stdev):
    """Return a copy of the structure with normally distributed random displacements applied to all positions.

    :param structure: the ``StructureData`` instance to rattle.
    :param stdev: ``Float`` with the standard deviation of the displacements in Angstrom.
    :returns: the rattled ``StructureData`` instance.
    """
    atoms = structure.get_ase()
    atoms.rattle(stdev=stdev.value)
    return orm.StructureData(ase=atoms)
