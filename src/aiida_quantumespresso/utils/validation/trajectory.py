# -*- coding: utf-8 -*-
"""Utilities for the validation of `TrajectoryData` content."""
from typing import List, Optional

from aiida.orm import TrajectoryData
import numpy


def verify_convergence_trajectory(
    trajectory: TrajectoryData,
    index: int = -1,
    threshold_forces: Optional[float] = None,
    threshold_stress: Optional[float] = None,
    reference_pressure: float = 0,
    fixed_coords: Optional[List[List[bool]]] = None
) -> bool:
    """Verify that the data of the given ``TrajectoryData`` is converged with respect to the given thresholds.

    Two are properties checked for convergence: forces and stress. If no threshold is specified for a property, the
    check is skipped. If all thresholds are successfully met, `True` is returned. If at least one of them fails, `False`
    is returned.

    The stress is converged if the absolute difference of its corresponding pressure and the reference pressure is
    within the given threshold.

    For calculations with fixed coordinates using the ``fixed_coords`` setting, the fixed coordinates can be provided.
    The forces that correspond to a fixed coordinate are ignored when comparing with the threshold.

    :param trajectory: the `TrajectoryData`
    :param index: the frame index of the trajectory data to check, default is `-1` meaning the last frame
    :param threshold_forces: the force threshold in Ry / bohr
    :param threshold_stress: the stress threshold in kbar
    :param reference_pressure: reference pressure in kbar
    :param fixed_coords: list of fixed coordinates in the calculation for each site
    :return: `True` if all thresholds are valid, `False` otherwise
    :raises ValueError: if any of the arrays or indices don't exist
    """
    if threshold_forces is not None:
        converged_forces = verify_convergence_forces(trajectory, index, threshold_forces, fixed_coords)
    else:
        converged_forces = True

    if threshold_stress is not None:
        converged_stress = verify_convergence_stress(trajectory, index, threshold_stress, reference_pressure)
    else:
        converged_stress = True

    return converged_forces and converged_stress


def verify_convergence_forces(
    trajectory: TrajectoryData,
    index: int = -1,
    threshold: Optional[float] = None,
    fixed_coords: Optional[List[List[bool]]] = None
) -> bool:
    """Verify that the `forces` of the given `TrajectoryData` are converged with respect to given threshold.

    For calculations with fixed coordinates using the ``fixed_coords`` setting, the fixed coordinates can be provided.
    The forces that correspond to a fixed coordinate are ignored when comparing with the threshold.

    :param trajectory: the `TrajectoryData`
    :param index: the frame index of the trajectory data to check, default is `-1` meaning the last frame
    :param threshold: the force threshold in Ry / bohr
    :param fixed_coords: list of fixed coordinates in the calculation for each site
    :return: `True` if threshold is valid, `False` otherwise
    :raises ValueError: if the `forces` array or given index does not exist
    """
    from qe_tools import CONSTANTS

    if threshold is None:
        return None

    threshold *= CONSTANTS.ry_to_ev / CONSTANTS.bohr_to_ang  # Convert to eV / Å

    try:
        abs_forces = numpy.abs(trajectory.get_array('forces')[index])
    except (KeyError, IndexError) as exception:
        raise ValueError('the `forces` array does not exist or the given index exceeds the length.') from exception

    if fixed_coords is not None:
        fixed_coords = numpy.array(fixed_coords)
        # Set the forces corresponding to fixed coordinates to zero, as they should not be checked versus threshold
        # Since `fixed_coords` is a list of lists of booleans, where `True` indicates that a coordinate is fixed, we can
        # invert this array and multiply it with the forces to set all fixed coordinates to zero.
        abs_forces *= numpy.invert(fixed_coords)

    return numpy.all(abs_forces < threshold)


def verify_convergence_stress(
    trajectory: TrajectoryData,
    index: int = -1,
    threshold: Optional[float] = None,
    reference_pressure: float = 0,
) -> bool:
    """Verify that the `stress` of the given `TrajectoryData` are converged with respect to given threshold.

    The stress is converged if the absolute difference of its corresponding pressure and the reference pressure is
    within the given threshold.

    :param trajectory: the `TrajectoryData`
    :param index: the frame index of the trajectory data to check, default is `-1` meaning the last frame
    :param threshold: the stress threshold in kbar
    :param reference_pressure: reference pressure in kbar
    :return: `True` if threshold is valid, `False` otherwise
    :raises ValueError: if the `stress` array or given index does not exist
    """
    if threshold is None:
        return None

    threshold /= 10.  # Convert to GPa
    reference_pressure /= 10.  # Convert to GPa

    try:
        stress = trajectory.get_array('stress')[index]
    except (KeyError, IndexError) as exception:
        raise ValueError('the `stress` array does not exist or the given index exceeds the length.') from exception

    pressure = numpy.trace(stress) / 3.

    return abs(pressure - reference_pressure) < threshold
