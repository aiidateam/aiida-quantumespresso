# -*- coding: utf-8 -*-
"""Utility class for handling the :py:meth:`data.hubbard_structure.HubbardStructureData`."""
from typing import List, Union
from pydantic import FilePath

from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

qe_translations = [
    [0, 0, 0],   [-1, -1, -1], [-1, -1, 0], [-1, -1, 1], [-1, 0, -1], [-1, 0, 0], [-1, 0, 1], 
    [-1, 1, -1], [-1, 1, 0],   [-1, 1, 1],  [0, -1, -1], [0, -1, 0],  [0, -1, 1], [0, 0, -1],
    [0, 0, 1],   [0, 1, -1],   [0, 1, 0],   [0, 1, 1],   [1, -1, -1], [1, -1, 0], [1, -1, 1],
    [1, 0, -1],  [1, 0, 0],    [1, 0, 1],   [1, 1, -1],  [1, 1, 0],   [1, 1, 1]
 ]
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
            self._hs = hubbard_structure
        else:
            raise ValueError('input is not of type `HubbardStructureData')

    @property
    def hubbard_structure(self) -> HubbardStructureData:
        """Return the HubbardStructureData."""
        return self._hs

    def get_hubbard_card(self) -> str:
        """Return QuantumESPRESSO `HUBBARD` input card for `pw.x`."""
        hubbard = self.hubbard_structure.hubbard
        hubbard_parameters = hubbard.parameters
        sites = self.hubbard_structure.sites
        natoms = len(sites)

        lines = [f'HUBBARD\t{hubbard.projectors}\n']

        for hp in hubbard_parameters:
            atom_i = sites[hp.atom_index].kind_name
            atom_j = sites[hp.neighbour_index].kind_name
            index_i = hp.atom_index +1 # QE indecis start from 1
            index_j = self._get_supercell_atomic_index(hp.neighbour_index, natoms, hp.translation) +1
            man_i = hp.atom_manifold
            man_j = hp.neighbour_manifold
            value = hp.value

            if hubbard.formulation not in ['dudarev', 'liechtenstein']:
                raise ValueError(f'Hubbard formulation {hubbard.formulation} is not implemented.')
            
            if hubbard.formulation=='liechtenstein':
                line = f'{pre}\t{atom_i}-{man_i} \t{value}'

            if hubbard.formulation=='dudarev':
                if hp.hubbard_type == 'J':
                    pre = 'J'
                elif atom_i == atom_j and hp.atom_manifold == hp.neighbour_manifold:
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

        .. note:: overrides any Hubbard stored information.
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
        if projectors.startswith(('(','[','{')):
            projectors = projectors[1:-1]

        # A sample of parsed array is: ['V', 'Co-3d', 'O-2p', '1', '4', '6.0']
        for data in hubbard_data:
            
            if data[0] == 'U':
                manifold = data[1].split('-')
                index = int(self.hubbard_structure._get_one_kind_index(manifold.pop(0))[0])
                manifold = '-'.join(manifold)
                args = (index, manifold, index, manifold, float(data[2]), [0,0,0], 'U')
            else:
                manifolds = []
                for i in [1,2]:
                    manifold = data[i].split('-')
                    manifold.pop(0) # removing atom name
                    manifolds.append('-'.join(manifold))

                index_i = int(data[3])
                index_j = int(data[4])
                found_i = False
                found_j = False

                for i, translation in enumerate(qe_translations):
                    if index_i - (i+1) * natoms <= 0 and not found_i:
                        index_0i = index_i - i * natoms - 1
                        found_i = True
                    if index_j - (i+1) * natoms <= 0 and not found_j:
                        index_0j = index_j - i * natoms - 1
                        translation_0 = translation
                        found_j = True

                args = (index_0i, manifolds[0], index_0j, manifolds[1], float(data[5]), translation_0, data[0])
            
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

        if not hubbard.formulation  == 'dudarev':
            raise ValueError('only `dudarev` formulation is implemented')

        card = '#\tAtom 1\tAtom 2\tHubbard V (eV)\n'

        for hp in hubbard_parameters:
            index_i = hp.atom_index +1 # QE indecis start from 1
            index_j = self._get_supercell_atomic_index(hp.neighbour_index, natoms, hp.translation) +1
            value = hp.value

            line = f'\t{index_i}\t{index_j}\t{value}'
            line += '\n'
            card += line

        return card

    def _get_supercell_atomic_index(self, index: int, num_sites: int, translation: List[Union[int,int,int]]) -> int:
        """Return the atomic index in 3x3x3 supercell.
        
        :param index: atomic index in unit cell
        :param num_sites: number of sites in structure
        :param translation: (3,) shape list of int referring to the translated atom in the 3x3x3 supercell
        
        :returns: atomic index in supercell standardized with the QuantumESPRESSO loop
        """
        return index + qe_translations.index(translation)*num_sites