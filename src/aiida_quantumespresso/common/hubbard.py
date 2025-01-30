# -*- coding: utf-8 -*-
"""Utility class and functions for HubbardStructureData."""
# pylint: disable=no-name-in-module, invalid-name
from typing import List, Literal, Tuple

from pydantic import BaseModel, conint, constr, field_validator

__all__ = ('HubbardParameters', 'Hubbard')


class HubbardParameters(BaseModel):
    """Class for describing onsite and intersite Hubbard interaction parameters.

    .. note: allowed manifold formats are:
            * {N}{L} (2 characters)
            * {N1}{L1}-{N2}{L2} (5 characters)

        N = quantum number (1,2,3,...); L = orbital letter (s,p,d,f,g,h)
    """

    atom_index: conint(strict=True, ge=0)
    """Atom index in the abstract structure."""

    atom_manifold: constr(strip_whitespace=True, to_lower=True, min_length=2, max_length=5)
    """Atom manifold (syntax is `3d`, `3d-2p`)."""

    neighbour_index: conint(strict=True, ge=0)
    """Neighbour index in the abstract structure."""

    neighbour_manifold: constr(strip_whitespace=True, to_lower=True, min_length=2, max_length=5)
    """Atom manifold (syntax is `3d`, `3d-2p`)."""

    translation: Tuple[conint(strict=True), conint(strict=True), conint(strict=True)]
    """Translation vector referring to the neighbour atom, (3,) shape list of ints."""

    value: float
    """Value of the Hubbard parameter, expessed in eV."""

    hubbard_type: Literal['Ueff', 'U', 'V', 'J', 'B', 'E2', 'E3']
    """Type of the Hubbard parameters used (`Ueff`, `U`, `V`, `J`, `B`, `E2`, `E3`)."""

    @field_validator('atom_manifold', 'neighbour_manifold')  # cls is mandatory to use
    def check_manifolds(cls, value):  # pylint: disable=no-self-argument, no-self-use
        """Check the validity of the manifold input.

        Allowed formats are:
            * {N}{L} (2 characters)
            * {N1}{L1}-{N2}{L2} (5 characters)

        N = quantum number (1,2,3,...); L = orbital letter (s,p,d,f,g,h)
        """
        length = len(value)
        if length not in [2, 5]:
            raise ValueError(f'invalid length ``{length}``. Only 2 or 5.')
        if length == 2:
            if not value[0] in [str(_ + 1) for _ in range(6)]:
                raise ValueError(f'invalid quantum number {value[0]}')
            if not value[1] in ['s', 'p', 'd', 'f', 'h']:
                raise ValueError(f'invalid manifold symbol {value[1]}')
        if length == 5:
            if not value[2] == '-':
                raise ValueError(f'the separator {value[0]} is not allowed. Only `-`')
            if not value[3] in [str(_ + 1) for _ in range(6)]:
                raise ValueError(f'the quantum number {value[0]} is not correct')
            if not value[4] in ['s', 'p', 'd', 'f', 'h']:
                raise ValueError(f'the manifold number {value[1]} is not correct')
        return value

    def to_tuple(self) -> Tuple[int, str, int, str, float, Tuple[int, int, int], str]:
        """Return the parameters as a tuple.

        The parameters have the following order:
            * atom_index
            * atom_manifold
            * neighbour_index
            * neighbour_manifold
            * value
            * translationr
            * hubbard_type
        """
        return (
            self.atom_index, self.atom_manifold, self.neighbour_index, self.neighbour_manifold, self.value,
            self.translation, self.hubbard_type
        )

    @staticmethod
    def from_tuple(hubbard_parameters: Tuple[int, str, int, str, float, Tuple[int, int, int], str]):
        """Return a ``HubbardParameters``  instance from a list.

        The parameters within the list must have the following order:
            * atom_index
            * atom_manifold
            * neighbour_index
            * neighbour_manifold
            * value
            * translation
            * hubbard_type
        """
        keys = [
            'atom_index',
            'atom_manifold',
            'neighbour_index',
            'neighbour_manifold',
            'value',
            'translation',
            'hubbard_type',
        ]
        return HubbardParameters(**dict(zip(keys, hubbard_parameters)))


class Hubbard(BaseModel):
    """Class for complete description of Hubbard interactions."""

    parameters: List[HubbardParameters]
    """List of :class:`~aiida_quantumespresso.common.hubbard.HubbardParameters`."""

    projectors: Literal['atomic',
                        'ortho-atomic',
                        'norm-atomic',
                        'wannier-functions',
                        'pseudo-potentials',
                        ] = 'ortho-atomic'
    """Name of the projectors used. Allowed values are:
        'atomic', 'ortho-atomic', 'norm-atomic', 'wannier-functions', 'pseudo-potentials'."""

    formulation: Literal['dudarev', 'liechtenstein'] = 'dudarev'
    """Hubbard formulation used. Allowed values are: 'dudarev', `liechtenstein`."""

    def to_list(self) -> List[Tuple[int, str, int, str, float, Tuple[int, int, int], str]]:
        """Return the Hubbard `parameters` as a list of lists.

        The parameters have the following order within each list:
            * atom_index
            * atom_manifold
            * neighbour_index
            * neighbour_manifold
            * value
            * translation
            * hubbard_type
        """
        return [params.to_tuple() for params in self.parameters]

    @staticmethod
    def from_list(
        parameters: List[Tuple[int, str, int, str, float, Tuple[int, int, int], str]],
        projectors: str = 'ortho-atomic',
        formulation: str = 'dudarev',
    ):
        """Return a :meth:`~aiida_quantumespresso.common.hubbard.Hubbard` instance from a list of tuples.

        Each list must contain the hubbard parameters in the following order:
            * atom_index
            * atom_manifold
            * neighbour_index
            * neighbour_manifold
            * value
            * translation
            * hubbard_type
        """
        parameters = [HubbardParameters.from_tuple(value) for value in parameters]
        return Hubbard(parameters=parameters, projectors=projectors, formulation=formulation)
