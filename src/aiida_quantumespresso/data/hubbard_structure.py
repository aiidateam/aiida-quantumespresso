# -*- coding: utf-8 -*-
"""Data plugin that represents a crystal structure with Hubbard parameters."""
import json
from typing import List, Tuple, Union

from aiida.orm import StructureData
import numpy as np
from pymatgen.core import Lattice, PeriodicSite

from aiida_quantumespresso.common.hubbard import Hubbard, HubbardParameters

__all__ = ('HubbardStructureData',)


class HubbardStructureData(StructureData):
    """Structure data containing code agnostic info on Hubbard parameters."""

    _hubbard_filename = 'hubbard.json'

    def __init__(
        self,
        cell: List[List[float]],
        sites: List[Tuple[str, str, Tuple[float, float, float]]],
        pbc: Tuple[bool, bool, bool] = (True, True, True),
        hubbard: Hubbard = None,
        **kwargs,
    ):
        """Set a ``HubbardStructureData`` instance.

        :param cell: (3,3) shape list of floats
        :param pbc: (3,) shape list in bools
        :param sites: list of lists, each of the shape [symbol, name, position],
            where position is a (3,) shape list of floats
        :param hubbard: a :class:`~aiida_quantumespresso.common.hubbard.Hubbard` istance
        """
        super().__init__(cell=cell, **kwargs)
        self.sites = sites
        self.pbc = pbc
        self.hubbard = Hubbard(parameters=[]) if hubbard is None else hubbard

    @property
    def sites(self):
        """Return the :attr:`~aiida.orm.nodes.data.structure.StructureData.sites`."""
        return super().sites

    @sites.setter
    def sites(self, values: List[Tuple[str, str, Tuple[float, float, float]]]):
        """Set the :attr:`~aiida.orm.nodes.data.structure.StructureData.sites`.

        :param values: list of sites, each as [symbol, name, (3,) shape list of positions]
        """
        self.clear_sites()
        for simbols, kind_name, position in values:
            self.append_atom(symbols=simbols, name=kind_name, position=position)

    @property
    def hubbard(self) -> Hubbard:
        """Get the `Hubbard` instance.

        :returns: a :class:`~aiida_quantumespresso.common.hubbard.Hubbard` instance.
        """
        # pylint: disable=not-context-manager
        with self.base.repository.open(self._hubbard_filename, mode='rb') as handle:
            return Hubbard.model_validate_json(json.load(handle))

    @hubbard.setter
    def hubbard(self, hubbard: Hubbard):
        """Set the full Hubbard information."""
        if not isinstance(hubbard, Hubbard):
            raise ValueError('the input is not of type `Hubbard`')

        serialized = json.dumps(hubbard.model_dump_json())
        self.base.repository.put_object_from_bytes(serialized.encode('utf-8'), self._hubbard_filename)

    @staticmethod
    def from_structure(
        structure: StructureData,
        hubbard: Union[Hubbard, None] = None,
    ):
        """Return an instance of ``HubbardStructureData`` from a ``StructureData`` node.

        :param structure: :class:`aiida.orm.StructureData` instance
        :param hubbard: :class:`~aiida_quantumespresso.common.hubbard.Hubbard` instance
        :returns: ``HubbardStructureData`` instance
        """
        sites = [[structure.get_kind(site.kind_name).symbol, site.kind_name, site.position] for site in structure.sites]
        cell = structure.cell
        pbc = structure.pbc

        return HubbardStructureData(cell=cell, pbc=pbc, sites=sites, hubbard=hubbard)

    def append_hubbard_parameter(
        self,
        atom_index: int,
        atom_manifold: str,
        neighbour_index: int,
        neighbour_manifold: str,
        value: float,
        translation: Tuple[int, int, int] = None,
        hubbard_type: str = 'Ueff',
    ):
        """Append a :class:`~aiida_quantumespresso.common.hubbard.HubbardParameters`.

        :param atom_index: atom index in unitcell
        :param atom_manifold: atomic manifold (e.g. 3d, 3d-2p)
        :param neighbour_index: neighbouring atom index in unitcell
        :param neighbour_manifold: neighbour manifold (e.g. 3d, 3d-2p)
        :param value: value of the Hubbard parameter, in eV
        :param translation: (3,) list of ints, describing the translation vector
            associated with the neighbour atom, defaults to None
        :param hubbard_type: hubbard type (U, V, J, ...), defaults to 'Ueff'
            (see :class:`~aiida_quantumespresso.common.hubbard.Hubbard` for full allowed values)
        """
        sites = [
            PeriodicSite(
                species=site.species,
                coords=site.coords,
                lattice=Lattice(self.cell, pbc=self.pbc),
                coords_are_cartesian=True
            ) for site in self.get_pymatgen().sites
        ]

        if any((atom_index > len(sites) - 1, neighbour_index > len(sites) - 1)):
            raise ValueError(
                'atom_index and neighbour_index must be within the range of the number of sites in the structure'
            )

        if translation is None:
            _, translation = sites[atom_index].distance_and_image(sites[neighbour_index])
            translation = np.array(translation, dtype=np.int64).tolist()

        hp_tuple = (atom_index, atom_manifold, neighbour_index, neighbour_manifold, value, translation, hubbard_type)
        parameters = HubbardParameters.from_tuple(hp_tuple)
        hubbard = self.hubbard

        if parameters in hubbard.parameters:
            hubbard.parameters.remove(parameters)

        hubbard.parameters.append(parameters)
        self.hubbard = hubbard

    def pop_hubbard_parameters(self, index: int):
        """Pop Hubbard parameters in the list.

        :param index: index of the Hubbard parameters to pop
        """
        hubbard = self.hubbard
        hubbard.parameters.pop(index)
        self.hubbard = hubbard

    def clear_hubbard_parameters(self):
        """Clear all the Hubbard parameters."""
        hubbard = self.hubbard
        hubbard.parameters = []
        self.hubbard = hubbard

    def initialize_intersites_hubbard(
        self,
        atom_name: str,
        atom_manifold: str,
        neighbour_name: str,
        neighbour_manifold: str,
        value: float = 1e-8,
        hubbard_type: str = 'V',
        use_kinds: bool = True,
    ):
        """Initialize and append intersite Hubbard values between an atom and its neighbour(s).

        .. note:: this only initialize the value between the first neighbour. In case
            `use_kinds` is False, all the possible combination of couples having
            kind  name equal to symbol are initialized.

        :param atom_name: atom name in unitcell
        :param atom_manifold: atomic manifold (e.g. 3d, 3d-2p)
        :param neighbour_index: neighbouring atom name in unitcell
        :param neighbour_manifold: neighbour manifold (e.g. 3d, 3d-2p)
        :param value: value of the Hubbard parameter, in eV
        :param hubbard_type: hubbard type (U, V, J, ...), defaults to 'V'
            (see :class:`~aiida_quantumespresso.common.hubbard.Hubbard` for full allowed values)
        :param use_kinds: whether to use kinds for initializing the parameters; when False, it
            initializes all the ``Kinds`` matching the ``atom_name``
        """
        sites = self.get_pymatgen_structure().sites

        function = self._get_one_kind_index if use_kinds else self._get_symbol_indices
        atom_indices = function(atom_name)
        neigh_indices = function(neighbour_name)

        if atom_indices is None or neigh_indices is None:
            raise ValueError('species or kind names not in structure')

        for atom_index in atom_indices:
            for neighbour_index in neigh_indices:
                _, translation = sites[atom_index].distance_and_image(sites[neighbour_index])
                translation = np.array(translation, dtype=np.int64).tolist()
                args = (
                    atom_index, atom_manifold, neighbour_index, neighbour_manifold, value, translation, hubbard_type
                )
                self.append_hubbard_parameter(*args)

    def initialize_onsites_hubbard(
        self,
        atom_name: str,
        atom_manifold: str,
        value: float = 1e-8,
        hubbard_type: str = 'Ueff',
        use_kinds: bool = True,
    ):
        """Initialize and append onsite Hubbard values of atoms with specific name.

        :param atom_name: atom name in unitcell
        :param atom_manifold: atomic manifold (e.g. 3d, 3d-2p)
        :param value: value of the Hubbard parameter, in eV
        :param hubbard_type: hubbard type (U, J, ...), defaults to 'Ueff'
            (see :class:`~aiida_quantumespresso.common.hubbard.Hubbard` for full allowed values)
        :param use_kinds: whether to use kinds for initializing the parameters; when False, it
            initializes all the ``Kinds`` matching the ``atom_name``
        """
        function = self._get_one_kind_index if use_kinds else self._get_symbol_indices
        atom_indices = function(atom_name)

        if atom_indices is None:
            raise ValueError('species or kind names not in structure')

        for atom_index in atom_indices:
            args = (atom_index, atom_manifold, atom_index, atom_manifold, value, [0, 0, 0], hubbard_type)
            self.append_hubbard_parameter(*args)

    def _get_one_kind_index(self, kind_name: str) -> List[int]:
        """Return the first site index matching with `kind_name`."""
        for i, site in enumerate(self.sites):
            if site.kind_name == kind_name:
                return [i]

    def _get_symbol_indices(self, symbol: str) -> List[int]:
        """Return one site index for each kind name matching symbol."""
        site_kindnames = self.get_site_kindnames()
        matching_kinds = [kind.name for kind in self.kinds if symbol in kind.symbol]

        return [site_kindnames.index(kind) for kind in matching_kinds]
