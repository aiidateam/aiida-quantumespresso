# -*- coding: utf-8 -*-
"""Utilities for `BandsData` nodes."""


def get_highest_occupied_band(bands, threshold=0.005):
    """Retun the index of the highest-occupied molecular orbital.

    The expected structure of the bands node is the following:

        * an array called `occupations`
        * with 3 dimensions if spin polarized, otherwise 2 dimensions
        * dimensions loop over, spin channel, kpoints, bands

    .. note::

        The threshold is used both as a limit below which a band is considered as unoccupied as well as a measure of
        numerical noise in the occupancies. The first band to have an occupancy below the threshold at all kpoints is
        marked as the LUMO. All subsequent bands are then expected to have an occupation that is less than twice the
        threshold. The reason for the factor two is that in the worst case when the LUMO has an occupation exactly equal
        to the threshold and a subsequent has twice that value, we can still consider that as not to break the
        requirement of all bands above the LUMO to be empty. It can be considered as having the exact same value (equal
        to the threshold) plus an numerical error also equal to the threshold.

    :param bands: the `BandsData` node
    :param threshold: raise a `ValueError` if the last band has an occupation above this threshold at any k-point
    :raises ValueError: if `bands` is not a `BandsData` node
    :raises ValueError: if `bands` does not contain the array `occupations`
    :raises ValueError: if `occupations` array has an invalid shape
    :raises ValueError: if any occupation above LUMO exceeds `2 * threshold`
    :raises ValueError: if the last band has an occupation above the threshold
    """
    from aiida.orm import BandsData
    from numpy import shape

    if not isinstance(bands, BandsData):
        raise ValueError(f'bands should be a `{BandsData.__name__}` node')

    try:
        occupations = bands.get_array('occupations')
    except KeyError as exception:
        raise ValueError('BandsData does not contain a `occupations` array') from exception

    lumo_indices = []

    # For spin-polarized calculations the `occupations` array should have 3 dimensions, otherwise just 2.
    if len(shape(occupations)) == 3:
        spin_channels = occupations
    elif len(shape(occupations)) == 2:
        spin_channels = [occupations]
    else:
        raise ValueError('invalid shape for `occupations` array')

    for l, spin_channel in enumerate(spin_channels):  # pylint: disable=invalid-name
        for k, kpoint in enumerate(spin_channel):

            lumo_index = None
            lumo_occupation = None

            for n, occupation in enumerate(kpoint):  # pylint: disable=invalid-name
                if lumo_index is not None:

                    # If the occupation of this band exceeds twice that of the threshold, it is considered to be not
                    # empty and since it comes after the LUMO, we raise
                    if occupation > 2 * threshold:
                        warning_args = [occupation, n, lumo_occupation, lumo_index, l, k]
                        raise ValueError('Occupation of {} at n={} after lumo lkn<{},{},{}>'.format(*warning_args))

                elif occupation < threshold:
                    lumo_index = n
                    lumo_occupation = occupation
                    lumo_indices.append(lumo_index)

            else:  # pylint: disable=useless-else-on-loop
                if kpoint[-1] >= threshold:
                    warning_args = [kpoint[-1], l, k, len(kpoint)]
                    raise ValueError('Occupation of {} at last band lkn<{},{},{}>'.format(*warning_args))

    # Note that the LUMO band indices are 0-indexed, so the actual band number is one higher, but the band number of the
    # HOMO is one lower than that, which therefore corresponds exactly to the 0-indexed LUMO index
    homo = max(lumo_indices)

    return homo
