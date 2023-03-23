# -*- coding: utf-8 -*-
"""Utility class for handling the :py:meth:`data.hubbard_structure.HubbardStructureData`."""
from math import floor
from typing import List, Union

from pydantic import FilePath  # pylint: disable=no-name-in-module

from aiida_quantumespresso.common.hubbard import Hubbard
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

qe_translations = [[0, 0, 0], [-1, -1, -1], [-1, -1, 0], [-1, -1, 1], [-1, 0, -1], [-1, 0, 0], [-1, 0, 1], [-1, 1, -1],
                   [-1, 1, 0], [-1, 1, 1], [0, -1, -1], [0, -1, 0], [0, -1, 1], [0, 0, -1], [0, 0, 1], [0, 1, -1],
                   [0, 1, 0], [0, 1, 1], [1, -1, -1], [1, -1, 0], [1, -1, 1], [1, 0, -1], [1, 0, 0], [1, 0, 1],
                   [1, 1, -1], [1, 1, 0], [1, 1, 1]]
# They can be generated with the following snippet:
#     from itertools import product
#     qe_translations = list(list(item) for item in product((-1, 0, 1), repeat=3))
#     first = qe_translations.pop(13)
#     qe_translations.insert(0, first)
#     return qe_translations


class HubbardUtils:
    """Utility class for handling `HubbardStructureData` for QuantumESPRESSO."""

    def __init__(
        self,
        hubbard_structure: HubbardStructureData,
    ):
        """Set a the `HubbardStructureData` to manipulate.

        :param hubbard_projectors: Hubbard projector type that are used for the occupations
        """
        if isinstance(hubbard_structure, HubbardStructureData):
            self._hubbard_structure = hubbard_structure
        else:
            raise ValueError('input is not of type `HubbardStructureData')

    @property
    def hubbard_structure(self) -> HubbardStructureData:
        """Return the HubbardStructureData."""
        return self._hubbard_structure

    def get_hubbard_card(self) -> str:
        """Return QuantumESPRESSO `HUBBARD` input card for `pw.x`."""
        hubbard = self.hubbard_structure.hubbard
        hubbard_parameters = hubbard.parameters
        sites = self.hubbard_structure.sites
        natoms = len(sites)

        lines = [f'HUBBARD\t{hubbard.projectors}\n']

        for param in hubbard_parameters:
            atom_i = sites[param.atom_index].kind_name
            atom_j = sites[param.neighbour_index].kind_name
            index_i = param.atom_index + 1  # QE indecis start from 1
            index_j = get_supercell_atomic_index(param.neighbour_index, natoms, param.translation) + 1
            man_i = param.atom_manifold
            man_j = param.neighbour_manifold
            value = param.value

            if hubbard.formulation not in ['dudarev', 'liechtenstein']:
                raise ValueError(f'Hubbard formulation {hubbard.formulation} is not implemented.')

            if hubbard.formulation == 'liechtenstein':
                line = f'{pre}\t{atom_i}-{man_i} \t{value}'

            is_intersite = is_intersite_hubbard(hubbard=hubbard)
            if hubbard.formulation == 'dudarev':
                if param.hubbard_type == 'J':
                    pre = 'J'
                elif is_intersite and atom_i == atom_j and param.atom_manifold == param.neighbour_manifold:
                    pre = 'U'
                else:
                    pre = 'V'

                if pre in ['U', 'J']:
                    line = f'{pre}\t{atom_i}-{man_i}\t{value}'
                else:
                    line = f'{pre}\t{atom_i}-{man_i}\t{atom_j}-{man_j}\t{index_i}\t{index_j}\t{value}'

            line += '\n'
            if line not in lines:
                lines.append(line)

        return ' '.join(lines)

    def print_hubbard_card(self):
        """Print the `HUBBARD` card."""
        print(self.get_hubbard_card())

    def parse_hubbard_dat(self, filepath: FilePath):
        """Parse the `HUBBARD.dat` of QuantumESPRESSO file associated to the current structure.

        This function is needed for parsing the HUBBARD.dat file generated in a `hp.x` calculation.

        .. note:: overrides current Hubbard information.
        """
        self.hubbard_structure.clear_hubbard_parameters()
        natoms = len(self.hubbard_structure.sites)
        hubbard_data = []

        with open(filepath, encoding='utf-8') as file:
            lines = file.readlines()
            for line in lines:
                if line.strip().split()[0] != '#':
                    hubbard_data.append(list(line.strip().split()))

        projectors = hubbard_data.pop(0)[1]
        if projectors.startswith(('(', '[', '{')):
            projectors = projectors[1:-1]

        # Samples of parsed array are:
        # ['U', 'Co-3d', '6.0']
        # ['V', 'Co-3d', 'O-2p', '1', '4', '6.0']
        for data in hubbard_data:

            if data[0] == 'U':
                manifold = data[1].split('-')
                index = int(self.hubbard_structure._get_one_kind_index(manifold.pop(0))[0])  # pylint: disable=protected-access
                manifold = '-'.join(manifold)
                args = (index, manifold, index, manifold, float(data[2]), [0, 0, 0], 'U')
            else:
                manifolds = []
                for i in [1, 2]:
                    manifold = data[i].split('-')
                    manifold.pop(0)  # removing atom name
                    manifolds.append('-'.join(manifold))

                # -1 because QE index starts from 1
                index_i, _ = get_supercell_site_properties(int(data[3]) - 1, natoms)  # pylint =
                index_j, tra = get_supercell_site_properties(int(data[4]) - 1, natoms)

                args = (index_i, manifolds[0], index_j, manifolds[1], float(data[5]), tra, data[0])

            self.hubbard_structure.append_hubbard_parameter(*args)

        hubbard = self.hubbard_structure.hubbard
        hubbard.projectors = projectors
        self.hubbard_structure.hubbard = hubbard

    def get_hubbard_file(self) -> str:
        """Return QuantumESPRESSO `parameters.in` data for `pw.x`."""
        hubbard = self.hubbard_structure.hubbard
        hubbard_parameters = hubbard.parameters
        sites = self.hubbard_structure.sites
        natoms = len(sites)

        if not hubbard.formulation == 'dudarev':
            raise ValueError('only `dudarev` formulation is implemented')

        card = '#\tAtom 1\tAtom 2\tHubbard V (eV)\n'

        for param in hubbard_parameters:
            index_i = param.atom_index + 1  # QE indecis start from 1
            index_j = get_supercell_atomic_index(param.neighbour_index, natoms, param.translation) + 1
            value = param.value

            line = f'\t{index_i}\t{index_j}\t{value}'
            line += '\n'
            card += line

        return card

    def reorder_atoms(self):
        """Reorder the atoms with with the kinds in the right order necessary for an ``hp.x`` calculation.

        An ``HpCalculation`` which restarts from a completed ``PwCalculation``, requires that the all
        Hubbard atoms appear first in  the atomic positions card of the ``PwCalculation`` input file.
        This order is based on the order of the kinds in the structure.
        So a suitable structure has all Hubbard kinds in the begining of kinds list.

        .. note:: overrides current ``HubbardStructureData``
        """
        structure = self.hubbard_structure  # current
        reordered = structure.clone()  # to be set at the end
        reordered.clear_kinds()

        hubbard = structure.hubbard.copy()
        parameters = hubbard.to_list()

        sites = structure.sites
        indices = get_hubbard_indices(hubbard=hubbard)
        hubbard_kinds = list(set([sites[index].kind_name for index in indices]))  # pylint: disable=consider-using-set-comprehension
        hubbard_kinds.sort(reverse=False)

        ordered_sites = []
        index_map = [index for index, site in enumerate(sites) if site.kind_name in hubbard_kinds]

        while hubbard_kinds:

            hubbard_kind = hubbard_kinds.pop()
            hubbard_sites = [s for s in sites if s.kind_name == hubbard_kind]
            remaining_sites = [s for s in sites if not s.kind_name == hubbard_kind]

            ordered_sites.extend(hubbard_sites)
            sites = remaining_sites

        # Extend the current site list with the remaining non-hubbard sites
        ordered_sites.extend(sites)

        for site in ordered_sites:

            if site.kind_name not in reordered.get_kind_names():
                kind = structure.get_kind(site.kind_name)
                reordered.append_kind(kind)

            reordered.append_site(site)

        reordered_parameters = []

        for parameter in parameters:
            new_parameter = parameter.copy()
            new_parameter[0] = index_map.index(parameter[0])
            new_parameter[2] = index_map.index(parameter[2])
            print(parameter, new_parameter)
            reordered_parameters.append(new_parameter)

        args = (reordered_parameters, hubbard.projectors, hubbard.formulation)
        hubbard = hubbard.from_list(*args)
        reordered.hubbard = hubbard

        self._hubbard_structure = reordered


def get_supercell_atomic_index(index: int, num_sites: int, translation: List[Union[int, int, int]]) -> int:
    """Return the atomic index in 3x3x3 supercell.

    :param index: atomic index in unit cell
    :param num_sites: number of sites in structure
    :param translation: (3,) shape list of int referring to the translated atom in the 3x3x3 supercell

    :returns: atomic index in supercell standardized with the QuantumESPRESSO loop
    """
    return index + qe_translations.index(translation) * num_sites


def get_supercell_site_properties(index: int, num_sites: int) -> List[Union[int, int, int]]:
    """Return the atomic index in unitcell and the associated translation from a 3x3x3 QuantumESPRESSO supercell index.

    :param index: atomic index
    :param num_sites: number of sites in structure
    """
    number = floor(index / num_sites)  # associated supercell number
    return (index - num_sites * number, qe_translations[number])


def get_hubbard_indices(hubbard: Hubbard) -> List[int]:
    """Return the set list of Hubbard indices."""
    atom_indecis = {parameters.atom_index for parameters in hubbard.parameters}
    neigh_indecis = {parameters.neighbour_index for parameters in hubbard.parameters}
    atom_indecis.update(neigh_indecis)
    return list(atom_indecis)


def is_intersite_hubbard(hubbard: Hubbard) -> bool:
    """Return whether `Hubbard` contains intersite interactions (+V)."""
    couples = [(param.atom_index != param.neighbour_index) for param in hubbard.parameters]
    return any(couples)
