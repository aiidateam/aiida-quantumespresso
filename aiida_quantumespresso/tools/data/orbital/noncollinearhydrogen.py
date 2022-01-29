# -*- coding: utf-8 -*-
"""A module defining hydrogen-like orbitals with non-collinear spin component."""

from aiida.common.exceptions import ValidationError
from aiida.tools.data.orbital.realhydrogen import RealhydrogenOrbital


def validate_spin(value):
    """Validate the value of the spin projection s_z."""
    if not isinstance(value, (int, float)):
        raise ValidationError('spin projection (s_z) must be float')

    return value


class NoncollinearHydrogenOrbital(RealhydrogenOrbital):
    """Orbitals for hydrogen, with non-collinear spin.

    The class largely follows the conventions used by wannier90 Object to
    handle the generation of real hydrogen orbitals and their hybrids, has
    methods for producing s, p, d, f, and sp, sp2, sp3, sp3d, sp3d2 hybrids.
    This method does not deal with the cannonical hydrogen orbitals which
    contain imaginary components. The orbitals described here are chiefly
    concerned with the symmetric aspects of the oribitals without the context
    of space. Therefore diffusitivity, position and atomic labels should be
    handled in the OrbitalData class. Following the notation of table 3.1, 3.2
    of Wannier90 user guide http://www.wannier.org/doc/user_guide.pdf A brief
    description of what is meant by each of these labels: :param radial_nodes:
    the number of radial nodes (or inflections) if no radial nodes are
    supplied, defaults to 0 :param angular_momentum: Angular quantum number,
    using real orbitals :param magnetic_number: Magnetic quantum number, using
    real orbitals :param spin: spin z-projection s_z The conventions regarding
    L and M correpsond to those used in wannier90 for all L greater than 0 the
    orbital is not hyrbridized see table 3.1 and for L less than 0 the orbital
    is hybridized see table 3.2. M then indexes all the possible orbitals from
    0 to 2L for L > 0 and from 0 to (-L) for L < 0.
    """

    _base_fields_required = RealhydrogenOrbital._base_fields_required

    _base_fields_optional = tuple(
        list(filter(lambda x: x[0] != 'spin', RealhydrogenOrbital._base_fields_optional)) + [
            ('spin', validate_spin, 0.0),
        ]
    )

    def __str__(self):
        """Printable representation of the orbital."""
        orb_dict = self.get_orbital_dict()
        try:
            orb_name = self.get_name_from_quantum_numbers(
                orb_dict['angular_momentum'], magnetic_number=orb_dict['magnetic_number']
            )
            pos = orb_dict['position']
            pos_string = f'{pos[0]:.4f},{pos[1]:.4f},{pos[2]:.4f}'
            orb = f"for kind {orb_dict['kind_name']}" if orb_dict['kind_name'] else ''
            out_string = f"r{orb_dict['radial_nodes']} {orb_name} (s_z={orb_dict['spin']}) orbital {orb} @ {pos_string}"
        except KeyError:
            # Should not happen, but we want it not to crash in __str__
            out_string = '(not all parameters properly set)'
        return out_string
