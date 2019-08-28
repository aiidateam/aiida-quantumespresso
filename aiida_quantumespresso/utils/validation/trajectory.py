# -*- coding: utf-8 -*-
"""Utilities for the validation of `TrajectoryData` content."""
from __future__ import absolute_import
from __future__ import division

import numpy


def verify_convergence_trajectory(trajectory, index=-1, threshold_forces=None, threshold_stress=None):
    """Verify that the data of the given `TrajectoryData` is converged with respect to the given thresholds.

    There are up to three properties controlled for convergence: total energy, forces and stress.
    If no threshold is specified the check is skipped. If any of the checks fails, either because the corresponding
    array, or index within that array, does not exist in the `TrajectoryData`, None is returned. If all thresholds are
    successfully met, `True` is returned, if at least one of them fails, `False` is returned.

    :param trajectory: the `TrajectoryData`
    :param index: the frame index of the trajectory data to check, default is `-1` meaning the last frame
    :param threshold_forces: the force threshold in Ry / bohr
    :param threshold_stress: the stress threshold in kbar
    :return: `None` if any of arrays or indices don't exist, `True` if all thresholds are valid, `False` otherwise
    """
    if threshold_forces is not None:
        converged_forces = verify_convergence_forces(trajectory, index, threshold_forces)
    else:
        converged_forces = True

    if threshold_stress is not None:
        converged_stress = verify_convergence_stress(trajectory, index, threshold_stress)
    else:
        converged_stress = True

    if converged_forces is None or converged_stress is None:
        return None

    return converged_forces & converged_stress


def verify_convergence_forces(trajectory, index=-1, threshold=None):
    """Verify that the `forces` of the given `TrajectoryData` are converged with respect to given threshold.

    :param trajectory: the `TrajectoryData`
    :param index: the frame index of the trajectory data to check, default is `-1` meaning the last frame
    :param threshold: the force threshold in Ry / bohr
    :return: `None` if the `forces` array or given index does not exist, `True` if threshold is valid, `False` otherwise
    """
    from qe_tools.constants import ry_to_ev, bohr_to_ang

    if threshold is None:
        return None

    threshold *= ry_to_ev / bohr_to_ang  # Convert to eV / Å

    try:
        forces = trajectory.get_array('forces')[index]
    except (KeyError, IndexError):
        return None

    return numpy.max(abs(forces)) < threshold


def verify_convergence_stress(trajectory, index=-1, threshold=None):
    """Verify that the `stress` of the given `TrajectoryData` are converged with respect to given threshold.

    :param trajectory: the `TrajectoryData`
    :param index: the frame index of the trajectory data to check, default is `-1` meaning the last frame
    :param threshold: the stress threshold in kbar
    :return: `None` if the `stress` array or given index does not exist, `True` if threshold is valid, `False` otherwise
    """
    if threshold is None:
        return None

    threshold /= 10.  # Convert to GPa

    try:
        stress = trajectory.get_array('stress')[index]
    except (KeyError, IndexError):
        return None

    return (numpy.trace(stress) / 3.) < threshold
