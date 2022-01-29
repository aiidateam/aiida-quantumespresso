# -*- coding: utf-8 -*-
"""A module defining hydrogen-like orbitals.

The orbitals are defined in the simultaneous eigen-basis of H, J^2, L^2, S^2
and J_z operators, useful for projecting densities with spin-orbit coupling.
"""

from aiida.common.exceptions import ValidationError
from aiida.tools.data.orbital.orbital import Orbital


def validate_j(value):
    """Validate the value of the total angular momentum."""
    if not isinstance(value, (int, float)):
        raise ValidationError('total angular momentum (j) must be float')

    if not value % 0.5 == 0.0:
        raise ValidationError('total angular momentum (j) must be a half integer')

    if value < 0.5 or value > 3.5:
        raise ValidationError('total angular momentum (j) must be between 0.5 and 3.5')

    return value


def validate_l(value):
    """Validate the value of the angular momentum."""
    if not isinstance(value, int):
        raise ValidationError('angular momentum (l) must be integer')

    if value < 0 or value > 3:
        raise ValidationError('angular momentum (l) must be between 0 and 3')

    return value


def validate_mj(value):
    """Validate the value of the magnetic number."""
    if not isinstance(value, (int, float)):
        raise ValidationError('magnetic number (m_j) must be float')

    # without knowing j we cannot validate the value of m_j. We will call an additional function
    # in the validate function

    return value


def validate_kind_name(value):
    """Validate the value of the kind_name."""
    if value is not None and not isinstance(value, str):
        raise ValidationError('kind_name must be a string')

    return value


def validate_n(value):
    """Validate the value of the number of radial nodes."""
    if not isinstance(value, int):
        raise ValidationError('number of radial nodes (n) must be integer')

    if value < 0 or value > 2:
        raise ValidationError('number of radial nodes (n) must be between 0 and 2')

    return value


class SpinorbitHydrogenOrbital(Orbital):
    """Orbitals for hydrogen with spin-orbit interaction.

    The orbital is defined in the common basis of H, J^2, L^2, S^2, J_z
    operators, indexed by the quantum numbers n, j, l, j_z.

    A brief description of what is meant by each of these labels:
    :param radial_nodes: the number of radial nodes (or inflections) if no radial
    nodes are supplied, defaults to 0
    :param angular_momentum: Angular quantum number l
    :param total_angular_momentum: Total angular momentum number j
    :param magnetic_number: Magnetic quantum number m_j

    The total angular momentum quantum number j takes values `|l-s| <= j <=
    l+s` in integer steps and the magnetic number takes values from -j to +j in
    integer steps.
    """

    _base_fields_required = tuple(
        list(Orbital._base_fields_required) +
        [('total_angular_momentum',
          validate_j), ('angular_momentum', validate_l), ('magnetic_number', validate_mj), ('radial_nodes', validate_n)]
    )

    _base_fields_optional = tuple(list(Orbital._base_fields_optional) + [
        ('kind_name', validate_kind_name, None),
    ])

    def __str__(self):
        """Printable representation of the orbital."""
        orb = self.get_orbital_dict()
        try:
            orb_name = f"j={orb['total_angular_momentum']},l={orb['angular_momentum']},m_j={orb['magnetic_number']}"
            position_string = f"{orb['position'][0]:.4f},{orb['position'][1]:.4f},{orb['position'][2]:.4f}"
            orbital = f"for kind {orb['kind_name']}" if orb['kind_name'] else ''
            out_string = f"r{orb['radial_nodes']} {orb_name} orbital {orbital} @ {position_string}"
        except KeyError:
            # Should not happen, but we want it not to crash in __str__
            out_string = '(not all parameters properly set)'
        return out_string

    def _validate_keys(self, input_dict):
        """Validate the keys otherwise raise ValidationError.

        Does basic validation from the parent followed by validations for the
        quantum numbers. Raises exceptions should the input_dict fail the
        valiation or if it contains any unsupported keywords. :param
        input_dict: the dictionary of keys to be validated :return
        validated_dict: a validated dictionary
        """
        validated_dict = super()._validate_keys(input_dict)

        # Validate m knowing the value of l
        total_angular_momentum = validated_dict['total_angular_momentum']  # j quantum number, must be there
        angular_momentum = validated_dict['angular_momentum']  # l quantum number, must be there
        accepted_range = [abs(angular_momentum - 0.5), angular_momentum + 0.5]
        if total_angular_momentum < min(accepted_range) or total_angular_momentum > max(accepted_range):
            raise ValidationError(
                f'the total angular momentum must be in the range [{min(accepted_range)}, {max(accepted_range)}]'
            )
        magnetic_number = validated_dict['magnetic_number']  # m quantum number, must be there
        accepted_range = [-total_angular_momentum, total_angular_momentum]
        if magnetic_number < min(accepted_range) or magnetic_number > max(accepted_range):
            raise ValidationError(
                f'the magnetic number must be in the range [{min(accepted_range)}, {max(accepted_range)}]'
            )

        return validated_dict
