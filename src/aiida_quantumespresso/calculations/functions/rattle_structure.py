"""Calculation function to apply random displacements to the atomic positions of a structure."""

import numpy as np
from aiida.engine import calcfunction


@calcfunction
def rattle_structure(structure, stdev):
    """Return a clone of the structure with normally distributed random displacements applied to all positions.

    The node is cloned rather than reconstructed, so any ``StructureData`` subclass (e.g.
    ``HubbardStructureData``) keeps its type, kind names and additional attributes.

    :param structure: the ``StructureData`` instance to rattle.
    :param stdev: ``Float`` with the standard deviation of the displacements in Angstrom.
    :returns: the rattled structure node.
    """
    rattled = structure.clone()
    positions = np.array([site.position for site in structure.sites])
    rng = np.random.RandomState(seed=42)  # same default seed as ``ase.Atoms.rattle``
    rattled.reset_sites_positions((positions + rng.normal(scale=stdev.value, size=positions.shape)).tolist())
    return rattled
