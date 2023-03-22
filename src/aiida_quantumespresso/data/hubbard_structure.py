# -*- coding: utf-8 -*-
"""Data plugin that represents a crystal structure with Hubbard parameters."""
import json
from typing import List, Union

from aiida.orm import StructureData
import numpy as np

from aiida_quantumespresso.common.hubbard import Hubbard, HubbardParameters


class HubbardStructureData(StructureData):
    """Structure data containing code independent info on Hubbard parameters."""

    _hubbard_filename = 'hubbard.json'

    def __init__(
        self,
        cell: List[List[float]],
        sites: List[Union[str, str, List[float]]],
        pbc: List[bool] = [True, True, True],
        hubbard: Hubbard = None,
        **kwargs,
    ):
        """Set a `HubbardStructureData` instance.

        :param cell: (3,3) shape list of ints
        :param pbc: (3,) shape list in bools
        :param sites: list of lists, each of the shape [symbol, name, position],
            where position is a (3,) shape list
        :param hubbard: a :py:meth:`~aiida_quantumespresso.common.hubbrd.Hubbard` istance
        """
        super().__init__(cell=cell, **kwargs)
        self.sites = sites
        self.pbc = pbc
        self.hubbard = Hubbard(parameters=[]) if hubbard is None else hubbard  # pylint: disable=dangerous-default-value

    @property
    def sites(self):
        """Return the `Sites`."""
        return super().sites

    @sites.setter
    def sites(self, values):
        """Set the `Sites`."""
        for simbols, kind_name, position in values:
            # Is it enough?
            self.append_atom(symbols=simbols, name=kind_name, position=position)

    @property
    def hubbard(self) -> Hubbard:
        """Get the `Hubbard` instance.

        :return: a :py:meth:`~aiida_quantumespresso.common.hubbard.Hubbard` instance.
        """
        with self.base.repository.open(self._hubbard_filename, mode='rb') as handle:
            return Hubbard.parse_raw(json.load(handle))

    @hubbard.setter
    def hubbard(self, hubbard: Hubbard):
        """Set the full Hubbard information."""
        if not isinstance(hubbard, Hubbard):
            raise ValueError('the input is not of type `Hubbard`')
        else:
            serialized = json.dumps(hubbard.json())
            self.base.repository.put_object_from_bytes(serialized.encode('utf-8'), self._hubbard_filename)

    @staticmethod
    def from_structure(
        structure: StructureData,
        hubbard: Hubbard = None,
    ):
        """Return an instance of HubbardStructureData from a StructureData object.

        :param structure: :meth:`aiida.orm.StructureData` instance
        :param parameters: the Hubbard parameters
        :param projectors: the Hubbard projector type
        :returns: HubbardStructureData instance
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
        translation: List[Union[int, int, int]] = None,
        hubbard_type: str = 'Ueff',
    ):
        """Append a Hubbard parameter."""
        pymat = self.get_pymatgen_structure()
        sites = pymat.sites

        if translation is None:
            _, translation = sites[atom_index].distance_and_image(sites[neighbour_index])
            translation = np.array(translation, dtype=np.int64).tolist()

        hp_list = [atom_index, atom_manifold, neighbour_index, neighbour_manifold, value, translation, hubbard_type]
        parameters = HubbardParameters.from_list(hp_list)
        hubbard = self.hubbard

        if parameters not in hubbard.parameters:
            hubbard.parameters.append(parameters)
            self.hubbard = hubbard

    def pop_hubbard_parameters(self, index: int):
        """Pop a Hubbard parameters."""
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
        hubbard_type: str = 'Ueff',
        use_kinds: bool = True,
    ):
        """Initialize and append intersite Hubbard values between an atom and its neighbour(s).

        .. note:: this only initialize the value between the first neighbour. In case
            `use_kinds` is False, all the possible combination of couples having
            kind  name equal to symbol are initialized.
        """
        sites = self.get_pymatgen_structure().sites

        function = self._get_one_kind_index if use_kinds else self._get_symbol_indecis
        atom_indecis = function(atom_name)
        neigh_indecis = function(neighbour_name)

        if atom_indecis is None or neigh_indecis is None:
            raise ValueError('species or kind names not in structure')

        for atom_index in atom_indecis:
            for neighbour_index in neigh_indecis:
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
        """Initialize and append onsite Hubbard values of atoms with specific name."""
        function = self._get_one_kind_index if use_kinds else self._get_symbol_indecis
        atom_indecis = function(atom_name)

        if atom_indecis is None:
            raise ValueError('species or kind names not in structure')

        for atom_index in atom_indecis:
            args = (atom_index, atom_manifold, atom_index, atom_manifold, value, [0, 0, 0], hubbard_type)
            self.append_hubbard_parameter(*args)

    def _get_one_kind_index(self, kind_name: str) -> List[int]:
        """Return the first site index matching with `kind_name`."""
        for i, site in enumerate(self.sites):
            if site.kind_name == kind_name:
                return [i]

    def _get_symbol_indecis(self, symbol: str) -> List[int]:
        """Return one site index for each kind name matching symbol."""
        site_kindnames = self.get_site_kindnames()
        matching_kinds = [kind.name for kind in self.kinds if symbol in kind.symbol]

        return [site_kindnames.index(kind) for kind in matching_kinds]
