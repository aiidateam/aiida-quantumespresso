# -*- coding: utf-8 -*-
"""Utilities for `BandsData` nodes."""
from __future__ import absolute_import


def get_highest_occupied_band(bands, threshold=0.01):
    """Retun the index of the highest-occupied molecular orbital.

    The expected structure of the bands node is the following:

    :param bands: the `BandsData` node
    :param threshold: raise a `ValueError` if the last band has an occupation above this threshold at any k-point
    :raises ValueError: if `bands` is not a `BandsData` node
    :raises ValueError: if `bands` does not contain the array `occupations`
    :raises ValueError: if `occupations` array has an invalid shape
    :raises ValueError: if the occupations are not monotonically decreasing at all k-points between LUMO and HOMO
    :raises ValueError: if the last band has an occupation above the threshold
    """
    from numpy import shape
    from aiida.orm import BandsData

    if not isinstance(bands, BandsData):
        raise ValueError('bands should be a `{}` node'.format(BandsData.__name__))

    try:
        occupations = bands.get_array('occupations')
    except KeyError:
        raise ValueError('BandsData does not contain a `occupations` array')

    lumo_indices = []

    # For spin-polarized calculations the `occupations` array should have 3 dimensions, otherwise just 2.
    if len(shape(occupations)) == 3:
        spin_channels = occupations
    elif len(shape(occupations)) == 2:
        spin_channels = [occupations]
    else:
        raise ValueError('invalid shape for `occupations` array')

    for l, spin_channel in enumerate(spin_channels):
        for k, kpoint in enumerate(spin_channel):

            lumo_index = None
            lumo_occupation = None

            for n, occupation in enumerate(kpoint):
                if lumo_index is not None:
                    if occupation > lumo_occupation:
                        warning_args = [occupation, n, lumo_occupation, lumo_index, l, k]
                        raise ValueError('Occupation of {} at n={} after lumo lkn<{},{},{}>'.format(*warning_args))
                elif occupation < threshold:
                    lumo_index = n
                    lumo_occupation = occupation
                    lumo_indices.append(lumo_index)
            else:
                if kpoint[-1] >= threshold:
                    warning_args = [kpoint[-1], l, k, len(kpoint)]
                    raise ValueError('Occupation of {} at last band lkn<{},{},{}>'.format(*warning_args))

    # Note that the LUMO band indices are 0-indexed, so the actual band number is one higher, but the band number of the
    # HOMO is one lower than that, which therefore corresponds exactly to the 0-indexed LUMO index
    homo = max(lumo_indices)

    return homo
