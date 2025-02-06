# -*- coding: utf-8 -*-
"""CalcFunction to compute the powder spectrum of a set of XANES spectra from ``XspectraCalculation``."""
import warnings

from aiida.common import ValidationError
from aiida.engine import calcfunction
from aiida.orm import XyData

warnings.warn(
    'This module is deprecated and will be removed soon as part of migrating XAS and XPS workflows to a new repository.'
    '\nThe new repository can be found at: https://github.com/aiidaplugins/aiida-qe-xspec.', FutureWarning
)


@calcfunction
def get_powder_spectrum(**kwargs):  # pylint: disable=too-many-statements
    """Combine the given spectra into a single "Powder" spectrum, representing the XAS of a powder sample.

    The function expects between 1 and 3 XyData nodes from ``XspectraCalculation`` whose
    polarisation vectors are the basis vectors of the original crystal structure (100, 010, 001).
    """
    spectra = [node for node in kwargs.values() if isinstance(node, XyData)]
    vectors = [node.creator.res['xepsilon'] for node in spectra]

    if len(vectors) > 3:
        raise ValidationError(f'Expected between 1 and 3 XyData nodes as input, but {len(spectra)} were given.')

    # If the system is isochoric (e.g. a cubic system) then the three basis vectors are
    # equal to each other, thus we simply return the
    if len(vectors) == 1:
        vector = vectors[0]
        vector_string = f'{float(vector[0])} {float(vector[1])} {float(vector[2])}'
        if vector not in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]:
            raise ValidationError(
                f'Polarisation vector ({vector_string}) does not correspond to a crystal basis vector (100, 010, 001)'
            )

        powder_spectrum = spectra[0]
        powder_x = powder_spectrum.get_x()[1]
        powder_y = powder_spectrum.get_y()[0][1]

        powder_data = XyData()
        powder_data.set_x(powder_x, 'energy', 'eV')
        powder_data.set_y(powder_y, 'sigma', 'n/a')

    # if the system is dichoric (e.g. a hexagonal system) then the A and B periodic
    # dimensions are equal to each other by symmetry, thus the powder spectrum is simply
    # the average of 2x the 1 0 0 eps vector and 1x the 0 0 1 eps vector
    if len(vectors) == 2:
        # Since the individual vectors are labelled, we can extract just the spectra needed
        # to produce the powder and leave the rest

        spectrum_a = None
        spectrum_c = None
        for vector, spectrum in zip(vectors, spectra):
            vector_string = f'{float(vector[0])} {float(vector[1])} {float(vector[2])}'
            if vector in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]:
                spectrum_a = spectrum
            elif vector == [0.0, 0.0, 1.0]:
                spectrum_c = spectrum
            else:
                raise ValidationError(
                    f'Polarisation vector ({vector_string}) does not correspond to a crystal basis vector'
                    ' (100, 010, 001)'
                )

        if spectrum_a and not spectrum_c:
            raise ValidationError(f'Found no polarisation vector for the C-axis ([0, 0, 1]), found instead: {vectors}.')
        powder_x = spectrum_a.get_x()[1]
        yvals_a = spectrum_a.get_y()[0][1]
        yvals_c = spectrum_c.get_y()[0][1]

        powder_y = ((yvals_a * 2) + yvals_c) / 3
        powder_data = XyData()
        powder_data.set_x(powder_x, 'energy', 'eV')
        powder_data.set_y(powder_y, 'sigma', 'n/a')

    # if the system is trichoric (e.g. a monoclinic system) then no periodic dimensions
    # are equal by symmetry, thus the powder spectrum is the average of the three basis
    # dipole vectors (1.0 0.0 0.0, 0.0 1.0 0.0, 0.0 0.0 1.0)
    if len(vectors) == 3:
        # Since the individual vectors are labelled, we can extract just the spectra needed to
        # produce the powder and leave the rest
        for vector, spectra in zip(vectors, spectra):
            vector_string = f'{float(vector[0])} {float(vector[1])} {float(vector[2])}'
            if vector == [1.0, 0.0, 0.0]:
                spectrum_a = spectra
            elif vector == [0.0, 1.0, 0.0]:
                spectrum_b = spectra
            elif vector == [0.0, 0.0, 1.0]:
                spectrum_c = spectra
            else:
                raise ValidationError(
                    f'Polarisation vector ({vector_string}) does not correspond to a crystal basis vector'
                    ' (100, 010, 001)'
                )

        powder_x = spectrum_a.get_x()[1]
        yvals_a = spectrum_a.get_y()[0][1]
        yvals_b = spectrum_b.get_y()[0][1]
        yvals_c = spectrum_c.get_y()[0][1]

        powder_y = (yvals_a + yvals_b + yvals_c) / 3

        powder_data = XyData()
        powder_data.set_x(powder_x, 'energy', 'eV')
        powder_data.set_y(powder_y, 'sigma', 'n/a')

    return powder_data
