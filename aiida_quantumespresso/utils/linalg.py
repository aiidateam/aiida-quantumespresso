# -*- coding: utf-8 -*-
"""Utilities related to linear algebra operations."""
import numpy


def are_matrices_equal(matrix_a, matrix_b, swap_sign_matrix_b=False, tolerance=1E-5):
    """Return whether matrix_a and matrix_b can be considered equal.

    The matrices will be cast to a numpy array and the difference will then be defined as the
    sum of all elements of the absolute difference between the two. If that sum is less than
    the tolerance, the function will return True. In all other cases it will be False

    :param matrix_a: a (nested) list
    :param matrix_b: a (nested) list
    :param swap_sign_matrix_b: when True will swap the sign of matrix_b before comparison
    :param tolerance: the tolerance within which the matrix difference is considered negligible
    :returns: True if the sum of the absolute difference elements is less than the tolerance, False otherwise
    """
    sign = -1 if swap_sign_matrix_b else 1

    return abs(numpy.array(matrix_a) - numpy.array(matrix_b) * sign).sum() < tolerance
