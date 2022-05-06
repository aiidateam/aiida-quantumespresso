# -*- coding: utf-8 -*-
"""Utilities for the validation of `TrajectoryData` content."""
import numpy


def verify_convergence_trajectory(
    trajectory, index=-1, threshold_forces=None, threshold_stress=None, reference_pressure=0.
):
    """Verify that the data of the given `TrajectoryData` is converged with respect to the given thresholds.

    There are up to three properties controlled for convergence: total energy, forces and stress.
    If no threshold is specified the check is skipped. If any of the checks fails, either because the corresponding
    array, or index within that array, does not exist in the `TrajectoryData`, None is returned. If all thresholds are
    successfully met, `True` is returned, if at least one of them fails, `False` is returned.

    The stress is converged if the absolute difference of its corresponding pressure and the reference pressure is
    within the given threshold.

    :param trajectory: the `TrajectoryData`
    :param index: the frame index of the trajectory data to check, default is `-1` meaning the last frame
    :param threshold_forces: the force threshold in Ry / bohr
    :param threshold_stress: the stress threshold in kbar
    :param reference_pressure: reference pressure in kbar
    :return: `True` if all thresholds are valid, `False` otherwise
    :raises ValueError: if any of the arrays or indices don't exist
    """
    if threshold_forces is not None:
        converged_forces = verify_convergence_forces(trajectory, index, threshold_forces)
    else:
        converged_forces = True

    if threshold_stress is not None:
        converged_stress = verify_convergence_stress(trajectory, index, threshold_stress, reference_pressure)
    else:
        converged_stress = True

    return converged_forces and converged_stress


def verify_convergence_forces(trajectory, index=-1, threshold=None):
    """Verify that the `forces` of the given `TrajectoryData` are converged with respect to given threshold.

    :param trajectory: the `TrajectoryData`
    :param index: the frame index of the trajectory data to check, default is `-1` meaning the last frame
    :param threshold: the force threshold in Ry / bohr
    :return: `True` if threshold is valid, `False` otherwise
    :raises ValueError: if the `forces` array or given index does not exist
    """
    from qe_tools import CONSTANTS

    if threshold is None:
        return None

    threshold *= CONSTANTS.ry_to_ev / CONSTANTS.bohr_to_ang  # Convert to eV / Å

    try:
        forces = trajectory.get_array('forces')[index]
    except (KeyError, IndexError) as exception:
        raise ValueError('the `forces` array does not exist or the given index exceeds the length.') from exception

    return numpy.max(abs(forces)) < threshold


def verify_convergence_stress(trajectory, index=-1, threshold=None, reference_pressure=0.):
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

    pressure = (numpy.trace(stress) / 3.)

    return abs(pressure - reference_pressure) < threshold
