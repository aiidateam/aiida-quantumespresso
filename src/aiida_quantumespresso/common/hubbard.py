# -*- coding: utf-8 -*-
"""Utility class and functions for HubbardStructureData."""
# pylint: disable=no-name-in-module, invalid-name
from typing import List, Literal, Union

from pydantic import BaseModel, conint, conlist, constr, validator

formulation_names = Literal['dudarev', 'liechtenstein']
type_names = Literal['Ueff', 'U', 'V', 'J', 'B', 'E2', 'E3']
projectors_name = Literal['atomic', 'ortho-atomic', 'norm-atomic', 'wannier-functions', 'pseudo-potentials',]


class HubbardParameters(BaseModel):
    """Class for describing onsite and intersite Hubbard interaction parameters."""

    atom_index: conint(strict=True, ge=0)
    atom_manifold: constr(strip_whitespace=True, to_lower=True, min_length=2, max_length=5)
    neighbour_index: conint(strict=True, ge=0)
    neighbour_manifold: constr(strip_whitespace=True, to_lower=True, min_length=2, max_length=5)
    translation: conlist(conint(strict=True), min_items=3, max_items=3)
    value: float
    hubbard_type: type_names  # default to 'U' ?

    @validator('atom_manifold', 'neighbour_manifold')  # cls is mandatory to use
    def check_manifolds(cls, value):  # pylint: disable=no-self-argument
        """Check the validity of the manifold input.

        Allowed formats are:
            * {N}{L} (2 characters)
            * {N1}{L1}-{N2}{L2} (5 characters)

        N = quantum number (1,2,3,...); L = orbital letter (s,p,d,f,g,h)
        """
        le = len(value)
        if le not in [2, 5]:
            raise ValueError(f'invalid length ``{le}``. Only 2 or 5.')
        if le == 2:
            if not value[0] in [str(_ + 1) for _ in range(6)]:
                raise ValueError(f'invalid quantum number {value[0]}')
            if not value[1] in ['s', 'p', 'd', 'f', 'h']:
                raise ValueError(f'invalid manifold symbol {value[1]}')
        if le == 5:
            if not value[2] == '-':
                raise ValueError(f'the separator {value[0]} is not allowed. Only `-`')
            if not value[3] in [str(_ + 1) for _ in range(6)]:
                raise ValueError(f'the quantum number {value[0]} is not correct')
            if not value[4] in ['s', 'p', 'd', 'f', 'h']:
                raise ValueError(f'the manifold number {value[1]} is not correct')
        return value

    def to_list(self):
        """Return the parameters as a list.

        The parameters have the following order:
            * atom_index
            * atom_manifold
            * neighbour_index
            * neighbour_manifold
            * value
            * translationr
            * hubbard_type
        """
        return [
            self.atom_index, self.atom_manifold, self.neighbour_index, self.neighbour_manifold, self.value,
            self.translation, self.hubbard_type
        ]

    @staticmethod
    def from_list(hubbard_parameters: list):
        """Return a `HubbardParameters` instance from a list.

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
        return HubbardParameters(dict(zip(keys, hubbard_parameters)))


class Hubbard(BaseModel):
    """Class for complete description of Hubbard interactions."""

    parameters: List[HubbardParameters]
    projectors: projectors_name = 'ortho-atomic'
    formulation: formulation_names = 'dudarev'

    def to_list(self) -> List[List[Union[int, str, int, str, float, list, str]]]:
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
        return [hp.to_list() for hp in self.parameters]

    @staticmethod
    def from_list(
        parameters: List[List[Union[int, str, int, str, float, list, str]]],
        projectors: str = 'ortho-atomic',
        formulation: str = 'dudarev',
    ):
        """Return a `Hubbard` instance from a list of lists.

        Each list must contain the hubbard parameters in the following order:
            * atom_index
            * atom_manifold
            * neighbour_index
            * neighbour_manifold
            * value
            * translation
            * hubbard_type
        """
        parameters = [HubbardParameters.from_list(value) for value in parameters]
        return Hubbard(parameters=parameters, projectors=projectors, formulation=formulation)
