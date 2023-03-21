# -*- coding: utf-8 -*-
"""Utility class and functions for HubbardStructureData."""
from typing import List, Union

from pydantic import BaseModel, ValidationError, conint, constr, validator

allowed_types_left = ['Dudarev', 'Liechtenstein']
allowed_types_right = ['U', 'V', 'J', 'Ueff', 'B', 'E2', 'E3']
allowed_projectors = [
    'atomic',
    'ortho-atomic',
    'norm-atomic',
    'wannier-functions',
    'pseudo-potentials',
]


class HubbardParameters(BaseModel):  # Or `HubbardInteraction` ?
    """Class for describing onsite and intersite Hubbard interaction parameters."""

    atom_index: conint(strict=True)
    atom_manifold: constr(strip_whitespace=True, to_lower=True, min_length=2, max_length=5)
    neighbour_index: conint(strict=True)
    neighbour_manifold: constr(strip_whitespace=True, to_lower=True, min_length=2, max_length=5)
    translation_vector: List[conint(strict=True)]
    hubbard_value: float  # just `value`?
    hubbard_type: str  # just `type`?

    @validator('atom_index', 'neighbour_index')
    def check_index(cls, value):
        """Check the positiveness of the atomic indecis."""
        if value < 0:
            raise ValueError('non integer index')
        return value

    @validator('translation_vector')
    def check_translation_vector(cls, value):
        """Check the validity of the traslation vector."""
        if len(value) != 3:
            raise ValueError('only 3 integers')
        return value

    @validator('atom_manifold', 'neighbour_manifold')
    def check_manifolds(cls, value):
        """Check the validity of the manifold input.

        Allowed formats are:
            * {N}{L} (2 characters)
            * {N1}{L1}-{N2}{L2} (5 characters)

        N = quantum number (1,2,3,...); L = orbital letter (s,p,d,f,g,h)
        """
        le = len(value)
        if le not in [2, 5]:
            raise ValueError(f'invalid length. Only 2 or 5.')
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

    @validator('hubbard_type')
    def check_hubbard_type(cls, value):
        """Check the validity of the `hubbard_type` input."""
        message = (
            'Hubbard type must be `{left}-{right}`, where the allowed values are:\n'\
            '* left: '+', '.join(allowed_types_left)+'\n'
            '* right: '+', '.join(allowed_types_right)
        )
        normalized = '-'.join((word.capitalize()) for word in value.split('-'))
        left = normalized.split('-')[0]
        try:
            right = normalized.split('-')[1]
        except IndexError:
            raise ValidationError(message)
        if left not in allowed_types_left or right not in allowed_types_right:
            raise ValueError(message)

        return normalized

    def to_list(self):
        """Return the parameters as a list.

        The parameters have the following order:
            * atom_index
            * atom_manifold
            * neighbour_index
            * neighbour_manifold
            * hubbard_value
            * translation_vector
            * hubbard_type
        """
        return [
            self.atom_index, self.atom_manifold, self.neighbour_index, self.neighbour_manifold, self.hubbard_value,
            self.translation_vector, self.hubbard_type
        ]

    @classmethod
    def from_list(cls, hubbard_parameters: list):
        """Return a `HubbardParameters` instance from a list.

        The parameters within the list must have the following order:
            * atom_index
            * atom_manifold
            * neighbour_index
            * neighbour_manifold
            * hubbard_value
            * translation_vector
            * hubbard_type
        """
        keys = [
            'atom_index',
            'atom_manifold',
            'neighbour_index',
            'neighbour_manifold',
            'hubbard_value',
            'translation_vector',
            'hubbard_type',
        ]
        return cls(**{key: value for key, value in zip(keys, hubbard_parameters)})


class HubbardProjectors(BaseModel):
    """Class defining the Hubbard projectors."""

    hubbard_projectors: str = 'ortho-atomic'  # only `projectors`?

    @validator('hubbard_projectors')
    def check_projectors(cls, value):
        """Check the validity of the `hubbard_parameters` input."""
        normalized = value.lower()
        if normalized not in allowed_projectors:
            raise ValueError('Hubbard projectors not in the allowed list: ' + ', '.join(allowed_projectors))
        return normalized


class Hubbard(HubbardProjectors):
    """Class for complete description of Hubbard interactions."""

    hubbard_parameters: List[HubbardParameters]  # only `parameters`?

    def to_list(self) -> List[List[Union[int, str, int, str, float, list, str]]]:
        """Return the `hubbard_parameters` as a list of lists.

        The parameters have the following order within each list:
            * atom_index
            * atom_manifold
            * neighbour_index
            * neighbour_manifold
            * hubbard_value
            * translation_vector
            * hubbard_type
        """
        return [hp.to_list() for hp in self.hubbard_parameters]

    @classmethod
    def from_list(cls, hubbard_parameters: List[List[Union[int, str, int, str, list, float, str]]]):
        """Return a `Hubbard` instance from a list of lists.

        Each list must contain the hubbard parameters in the following order:
            * atom_index
            * atom_manifold
            * neighbour_index
            * neighbour_manifold
            * hubbard_value
            * translation_vector
            * hubbard_type

        .. note:: the `hubbard_projectors` cannot be specified directly from this method
        """
        return cls(hubbard_parameters=[HubbardParameters.from_list(value) for value in hubbard_parameters])
