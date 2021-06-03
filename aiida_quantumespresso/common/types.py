# -*- coding: utf-8 -*-
"""Module with common data types."""
import enum


class ElectronicType(enum.Enum):
    """Enumeration to indicate the electronic type of a system."""

    METAL = 'metal'
    INSULATOR = 'insulator'
    AUTOMATIC = 'automatic'


class RelaxType(enum.Enum):
    """Enumeration of known relax types."""

    NONE = 'none'  # All degrees of freedom are fixed, essentially performs single point SCF calculation
    POSITIONS = 'positions'  # Only the atomic positions are relaxed, cell is fixed
    VOLUME = 'volume'  # Only the cell volume is optimized, cell shape and atoms are fixed
    SHAPE = 'shape'  # Only the cell shape is optimized at a fixed volume and fixed atomic positions
    CELL = 'cell'  # Only the cell is optimized, both shape and volume, while atomic positions are fixed
    POSITIONS_VOLUME = 'positions_volume'  # Same as `VOLUME` but atomic positions are relaxed as well
    POSITIONS_SHAPE = 'positions_shape'  # Same as `SHAPE`  but atomic positions are relaxed as well
    POSITIONS_CELL = 'positions_cell'  # Same as `CELL`  but atomic positions are relaxed as well


class SpinType(enum.Enum):
    """Enumeration to indicate the spin polarization type of a system."""

    NONE = 'none'
    COLLINEAR = 'collinear'
    NON_COLLINEAR = 'non_collinear'
    SPIN_ORBIT = 'spin_orbit'
