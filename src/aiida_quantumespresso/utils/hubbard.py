# -*- coding: utf-8 -*-
"""Utility class for handling the :class:`aiida_quantumespresso.data.hubbard_structure.HubbardStructureData`."""
# pylint: disable=no-name-in-module
from itertools import product
import os
from typing import Dict, List, Tuple, Union

from aiida.orm import StructureData
import numpy as np

from aiida_quantumespresso.common.hubbard import Hubbard
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData

__all__ = (
    'HubbardUtils',
    'initialize_hubbard_parameters',
    'get_supercell_atomic_index',
    'get_index_and_translation',
    'get_hubbard_indices',
    'is_intersite_hubbard',
    'max_number_of_neighbours',
)

QE_TRANSLATIONS = list(list(item) for item in product((-1, 0, 1), repeat=3))
first = QE_TRANSLATIONS.pop(13)
QE_TRANSLATIONS.insert(0, first)
QE_TRANSLATIONS = tuple(tuple(item) for item in QE_TRANSLATIONS)


class HubbardUtils:
    """Utility class for handling `HubbardStructureData` for QuantumESPRESSO."""

    def __init__(
        self,
        hubbard_structure: HubbardStructureData,
    ):
        """Set a the `HubbardStructureData` to manipulate."""
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
            index_i = param.atom_index + 1  # QE indices start from 1
            index_j = get_supercell_atomic_index(param.neighbour_index, natoms, param.translation) + 1
            man_i = param.atom_manifold
            man_j = param.neighbour_manifold
            value = param.value

            if hubbard.formulation not in ['dudarev', 'liechtenstein']:
                raise ValueError(f'Hubbard formulation {hubbard.formulation} is not implemented.')

            if hubbard.formulation == 'liechtenstein':
                line = f'{pre}\t{atom_i}-{man_i} \t{value}'

            # This variable is to meet QE implementation. If intersite interactions
            # (+V) are present, onsite parameters might not be relabelled by the ``hp.x``
            # code, causing a subsequent ``pw.x`` calculation to crash. That is,
            # we need to avoid writing "U Co-3d 5.0", but instead "V Co-3d Co-3d 1 1 5.0".
            is_intersite = is_intersite_hubbard(hubbard=hubbard)
            if hubbard.formulation == 'dudarev':
                if param.hubbard_type == 'J':
                    pre = 'J'
                elif not is_intersite and atom_i == atom_j and param.atom_manifold == param.neighbour_manifold:
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

    def parse_hubbard_dat(self, filepath: Union[str, os.PathLike]):
        """Parse the `HUBBARD.dat` of QuantumESPRESSO file associated to the current structure.

        This function is needed for parsing the HUBBARD.dat file generated in a `hp.x` calculation.

        .. note:: overrides current Hubbard information.

        :param filepath: the filepath of the *HUBBARD.dat* to parse
        """
        self.hubbard_structure.clear_hubbard_parameters()
        natoms = len(self.hubbard_structure.sites)
        hubbard_data = []

        with open(filepath, encoding='utf-8') as handle:
            lines = handle.readlines()
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
                args = (index, manifold, index, manifold, float(data[2]), (0, 0, 0), 'U')
            else:
                manifolds = []
                for i in [1, 2]:
                    manifold = data[i].split('-')
                    manifold.pop(0)  # removing atom name
                    manifolds.append('-'.join(manifold))

                # -1 because QE index starts from 1
                index_i, _ = get_index_and_translation(int(data[3]) - 1, natoms)
                index_j, tra = get_index_and_translation(int(data[4]) - 1, natoms)

                args = (index_i, manifolds[0], index_j, manifolds[1], float(data[5]), tuple(tra), data[0])

            self.hubbard_structure.append_hubbard_parameter(*args)

        hubbard = self.hubbard_structure.hubbard
        hubbard.projectors = projectors
        self.hubbard_structure.hubbard = hubbard

    def get_hubbard_file(self) -> str:
        """Return QuantumESPRESSO ``parameters.in`` data for ``pw.x```."""
        hubbard = self.hubbard_structure.hubbard
        hubbard_parameters = hubbard.parameters
        sites = self.hubbard_structure.sites
        natoms = len(sites)

        if not hubbard.formulation == 'dudarev':
            raise ValueError('only `dudarev` formulation is implemented')

        card = '#\tAtom 1\tAtom 2\tHubbard V (eV)\n'

        for param in hubbard_parameters:
            index_i = param.atom_index + 1  # QE indices start from 1
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
        from copy import deepcopy

        structure = self.hubbard_structure  # current
        reordered = structure.clone()  # to be set at the end
        reordered.clear_kinds()

        hubbard = structure.hubbard.model_copy()
        parameters = hubbard.to_list()

        sites = structure.sites
        indices = get_hubbard_indices(hubbard=hubbard)
        hubbard_kinds = list(set(sites[index].kind_name for index in indices))
        hubbard_kinds.sort(reverse=False)

        ordered_sites = []

        # We define a map from ``index`` to ``site specifications``. We need the complete
        # specification, as we will loose track of the index ordering with the following shuffle.
        # The ``index`` are needed for re-indexing later the hubbard parameters.
        index_map = {index: site.get_raw() for index, site in enumerate(sites) if site.kind_name in hubbard_kinds}

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

        # Reordered site map, to match with raw sites of ``index_map``
        site_map = [site.get_raw() for site in reordered.sites]

        for parameter in parameters:
            new_parameter = list(deepcopy(parameter))
            new_parameter[0] = site_map.index(index_map[parameter[0]])  # atom index
            new_parameter[2] = site_map.index(index_map[parameter[2]])  # neighbour index
            reordered_parameters.append(new_parameter)

        # making sure we keep track of the other info as well
        args = (reordered_parameters, hubbard.projectors, hubbard.formulation)
        hubbard = hubbard.from_list(*args)
        reordered.hubbard = hubbard

        self._hubbard_structure = reordered

    def is_to_reorder(self) -> bool:
        """Return whether the atoms should be reordered for an ``hp.x`` calculation."""
        indices = get_hubbard_indices(self.hubbard_structure.hubbard)
        indices.sort()

        return indices != list(range(len(indices)))

    def get_hubbard_for_supercell(self, supercell: StructureData, thr: float = 1e-3) -> HubbardStructureData:
        """Return the ``HubbbardStructureData`` for a supercell.

        .. note:: the two structure need to be commensurate (no rigid rotations)

        .. warning:: **always check** that the energy calculation of a pristine supercell
            structure obtained through this method is the same as the unitcell (within numerical noise)

        :returns: a new ``HubbbardStructureData`` with all the mapped Hubbard parameters
        """
        uc_pymat = self.hubbard_structure.get_pymatgen_structure()
        sc_pymat = supercell.get_pymatgen_structure()
        uc_positions = uc_pymat.cart_coords  # positions in Cartesian coordinates
        sc_positions = sc_pymat.cart_coords
        uc_cell = uc_pymat.lattice.matrix
        uc_cell_inv = np.linalg.inv(uc_cell)

        hubbard = self.hubbard_structure.hubbard
        sc_hubbard_parameters = []

        # Dumb, but fairly safe, way to map all the hubbard parameters.
        # The idea is to map for each interaction in unitcell the
        # correspective one in supercell matching all the positions.
        for param in hubbard.parameters:
            # i -> atom_index | j -> neighbour_index
            uc_i_position = uc_positions[param.atom_index]
            uc_j_position = uc_positions[param.neighbour_index]
            sc_i_indices, sc_j_index, sc_i_translations = [], [], []

            # Each atom in supercell is matched if a unitcell
            # translation vector is found.
            for i, position in enumerate(sc_positions):
                translation = np.dot(position - uc_i_position, uc_cell_inv)
                translation_int = np.rint(translation)
                if np.all(np.isclose(translation, translation_int, thr)):
                    sc_i_translations.append(translation_int.tolist())
                    sc_i_indices.append(i)

                translation = np.dot(position - uc_j_position, uc_cell_inv)
                translation_int = np.rint(translation)
                if np.all(np.isclose(translation, translation_int, thr)):
                    uc_j_translation = np.array(translation_int)
                    sc_j_index = i

            # The position of the neighbour must be still translated;
            # This might happen in the supercell itself, or outside, thus
            # we neeed to recompute its position and its translation vector in supercell.
            sc_j_position = sc_positions[sc_j_index]

            for sc_i_index, sc_i_translation in zip(sc_i_indices, sc_i_translations):
                j_position = sc_j_position + np.dot(sc_i_translation - uc_j_translation + param.translation, uc_cell)
                local_site = sc_pymat.get_sites_in_sphere(pt=j_position, r=thr)[0]  # pymatgen PeriodicSite

                sc_hubbard_parameter = [
                    int(sc_i_index),
                    param.atom_manifold,
                    int(local_site.index),  # otherwise the class validator complains
                    param.neighbour_manifold,
                    param.value,
                    np.array(local_site.image, dtype=np.int64).tolist(),  # otherwise the class validator complains
                    param.hubbard_type,
                ]

                sc_hubbard_parameters.append(sc_hubbard_parameter)

        args = (sc_hubbard_parameters, hubbard.projectors, hubbard.formulation)
        new_hubbard = Hubbard.from_list(*args)

        return HubbardStructureData.from_structure(structure=supercell, hubbard=new_hubbard)

    def get_interacting_pairs(self) -> Dict[str, List[str]]:
        """Return tuple of kind name interaction pairs.

        :returns: dictionary of onsite kinds with a list of V kinds
        """
        pairs_dict = dict()
        sites = self.hubbard_structure.sites
        parameters = self.hubbard_structure.hubbard.parameters

        onsite_indices = [p.atom_index for p in parameters if p.atom_index == p.neighbour_index]

        for index in onsite_indices:
            name = sites[index].kind_name
            if name not in pairs_dict:
                pairs_dict[name] = []

        for parameter in parameters:
            onsite_name = sites[parameter.atom_index].kind_name
            neigh_name = sites[parameter.neighbour_index].kind_name

            if onsite_name in pairs_dict and onsite_name != neigh_name:
                if not neigh_name in pairs_dict[onsite_name]:
                    pairs_dict[onsite_name].append(neigh_name)

        return pairs_dict

    def get_pairs_radius(
        self,
        onsite_index: int,
        neighbours_names: List[str],
        number_of_neighbours: int,
        radius_max: float = 7.0,
        thr: float = 1.0e-2,
    ) -> Tuple[float, float]:
        """Return the minimum and maximum radius of the first neighbours of the onsite site.

        :param onsite_index: index in the structure of the onsite Hubbard atom
        :param neighbours_names: kind names of the neighbours
        :param number_of_neighbours: number of neighbours coming to select
        :param radius_max: maximum radius (in Angstrom) to use for looking for neighbours
        :param thr: threshold (in Angstrom) for defining the shells
        :return: (radius min +thr, radius max -thr) defining the shells containing only the first neighbours
        """
        rmin = 0
        pymat = self.hubbard_structure.get_pymatgen_structure()
        _, neigh_indices, _, distances = pymat.get_neighbor_list(sites=[pymat[onsite_index]], r=radius_max)

        sort = np.argsort(distances)
        neigh_indices = neigh_indices[sort]
        distances = distances[sort]

        count = 0
        for i in range(len(neigh_indices)):  # pylint: disable=consider-using-enumerate
            index = i
            if self.hubbard_structure.sites[neigh_indices[i]].kind_name in neighbours_names:
                rmin = max(rmin, distances[i])
                count += 1
            if count == number_of_neighbours:
                break

        rmax = distances[index + 1]
        if np.abs(rmax - rmin) < thr:
            rmax = rmin + 2 * rmin

        return rmin + thr, rmax - thr

    def get_intersites_radius(
        self,
        nn_finder: str = 'crystal',
        nn_inputs: Union[Dict, None] = None,
        radius_max: float = 7.0,
        thr: float = 1.0e-2,
        **_,
    ) -> float:
        """Return the radius (in Angstrom) for intersites from nearest neighbour finders.

        It peforms a nearest neighbour analysis (via pymatgen modules) to find the first inersite
        neighbours for all the onsite atoms. A radius is returned which can be used to
        run an ``hp.x`` calculation. Such radius defines a shell including only the first
        neighbours of each onsite Hubbard atom.

        :param nn_finder: string defining the nearest neighbour finder; options are:
            * `crystal`: use :class:`pymatgen.analysis.local_env.CrystalNN`
            * `voronoi`: use :class:`pymatgen.analysis.local_env.VoronoiNN`
        :param nn_inputs: inputs for the nearest neighbours finder; when None, standard inputs
            are used to find geometric first neighbours (recommended)
        :param radius_max: max radius where to look for neighbouring atoms, in Angstrom
        :param thr: threshold (in Angstrom) for defining the shells
        :return: radius defining the shell containing only the first neighbours
        """
        import warnings

        from pymatgen.analysis.local_env import CrystalNN, VoronoiNN

        rmin, rmax = 0.0, radius_max

        if nn_finder not in ['crystal', 'voronoi']:
            raise ValueError('`nn_finder` must be either `cyrstal` or `voronoi`')
        if nn_inputs is None:
            if nn_finder == 'crystal':
                nn_inputs = {'distance_cutoffs': None, 'x_diff_weight': 0, 'porous_adjustment': False}
            if nn_finder == 'voronoi':
                nn_inputs = {'tol': 0.1, 'cutoff': radius_max}

        voronoi = CrystalNN(**nn_inputs) if nn_finder == 'crystal' else VoronoiNN(**nn_inputs)  # pylint: disable=unexpected-keyword-arg

        sites = self.hubbard_structure.sites
        name_to_specie = {kind.name: kind.symbol for kind in self.hubbard_structure.kinds}
        pymat = self.hubbard_structure.get_pymatgen_structure()
        pairs = self.get_interacting_pairs()

        for i, site in enumerate(sites):
            if site.kind_name in pairs:
                neigh_species = voronoi.get_cn_dict(pymat, i)  # e.g. {'O': 4, 'S': 2, ...}
                number_of_neighs = 0

                for neigh_name in pairs[site.kind_name]:
                    specie = name_to_specie[neigh_name]

                    if specie in neigh_species:
                        number_of_neighs += neigh_species[specie]
                        neigh_species.pop(specie)  # avoid 'duplicating' same specie but different (kind) name

                rmin_, rmax_ = self.get_pairs_radius(i, pairs[site.kind_name], number_of_neighs, radius_max, thr)

                rmin = max(rmin_, rmin)  # we want the largest to include them all
                rmax = min(rmax_, rmax)  # we want the smallest to check whether such radius exist

        if rmin > rmax:
            warnings.warn('A common radius seems to not exist! Try lowering `thr`.')

        return min(rmin, rmax)

    def get_intersites_list(
        self,
        nn_finder: str = 'crystal',
        nn_inputs: Union[Dict, None] = None,
        radius_max: float = 7.0,
        **_,
    ) -> List[Tuple[int, int, Tuple[int, int, int]]]:
        """Return the list of intersites from nearest neighbour finders.

        It peforms a nearest neighbour analysis (via pymatgen modules) to find the first inersite
        neighbours for all the onsite atoms. A list is returned with all the nearest neighbours
        providing all the information about the couples indices and the associated trasnaltion vector.
        Also on-site information is included.

        :param nn_finder: string defining the nearest neighbour finder; options are:
            * `crystal`: use :class:`pymatgen.analysis.local_env.CrystalNN`
            * `voronoi`: use :class:`pymatgen.analysis.local_env.VoronoiNN`
        :param nn_inputs: inputs for the nearest neighbours finder; when None, standard inputs
            are used to find geometric first neighbours (recommended)
        :param radius_max: max radius where to look for neighbouring atoms, in Angstrom
        :param thr: threshold (in Angstrom) for defining the shells
        :return: list of lists, each having (atom index, neighbouring index, translation vector)
        """
        from pymatgen.analysis.local_env import CrystalNN, VoronoiNN

        if nn_finder not in ['crystal', 'voronoi']:
            raise ValueError('`nn_finder` must be either `cyrstal` or `voronoi`')

        if nn_inputs is None:
            if nn_finder == 'crystal':
                nn_inputs = {'distance_cutoffs': None, 'x_diff_weight': 0, 'porous_adjustment': False}
            if nn_finder == 'voronoi':
                nn_inputs = {'tol': 0.1, 'cutoff': radius_max}

        voronoi = CrystalNN(**nn_inputs) if nn_finder == 'crystal' else VoronoiNN(**nn_inputs)  # pylint: disable=unexpected-keyword-arg

        sites = self.hubbard_structure.sites
        name_to_specie = {kind.name: kind.symbol for kind in self.hubbard_structure.kinds}
        pymat = self.hubbard_structure.get_pymatgen_structure()
        pairs = self.get_interacting_pairs()
        neigh_list = []

        for i, site in enumerate(sites):  # pylint: disable=too-many-nested-blocks
            if site.kind_name in pairs:
                neigh_species = voronoi.get_cn_dict(pymat, i)  # e.g. {'O': 4, 'S': 2, ...}
                neigh_list.append([i, i, (0, 0, 0)])

                for neigh_name in pairs[site.kind_name]:
                    try:  # we want an API that allows for specifying neighbours that are not present
                        specie = name_to_specie[neigh_name]

                        count = 0
                        if specie in neigh_species:
                            _, neigh_indices, images, distances = pymat.get_neighbor_list(
                                sites=[pymat[i]], r=radius_max
                            )

                            sort = np.argsort(distances)
                            neigh_indices = neigh_indices[sort]
                            images = images[sort]

                            for index, image in zip(neigh_indices, images):
                                if pymat[index].specie.name == specie:
                                    neigh_list.append([
                                        i,
                                        int(index),
                                        tuple(np.array(image, dtype=np.int64).tolist()),
                                    ])
                                    count += 1

                                if count >= neigh_species[specie]:
                                    break
                    except KeyError:
                        pass

        return neigh_list

    def get_max_number_of_neighbours(
        self,
        nn_finder: str = 'crystal',
        nn_inputs: Union[Dict, None] = None,
        radius_max: float = 7.0,
        **_,
    ) -> list:
        """Return the maximum number of nearest neighbours, aslo counting the non-interacting ones.

        It peforms a nearest neighbour analysis (via pymatgen modules) to find the first inersite
        neighbours for all the onsite atoms. A list is returned with all the nearest neighbours
        providing all the information about the couples indices and the associated trasnaltion vector.
        Also on-site information is included.

        :param nn_finder: string defining the nearest neighbour finder; options are:
            * `crystal`: use :class:`pymatgen.analysis.local_env.CrystalNN`
            * `voronoi`: use :class:`pymatgen.analysis.local_env.VoronoiNN`
        :param nn_inputs: inputs for the nearest neighbours finder; when None, standard inputs
            are used to find geometric first neighbours (recommended)
        :param radius_max: max radius where to look for neighbouring atoms, in Angstrom
        :param thr: threshold (in Angstrom) for defining the shells
        :return: list of lists, each having (atom index, neighbouring index, translation vector)
        """
        from pymatgen.analysis.local_env import CrystalNN, VoronoiNN

        if nn_finder not in ['crystal', 'voronoi']:
            raise ValueError('`nn_finder` must be either `cyrstal` or `voronoi`')

        if nn_inputs is None:
            if nn_finder == 'crystal':
                nn_inputs = {'distance_cutoffs': None, 'x_diff_weight': 0, 'porous_adjustment': False}
            if nn_finder == 'voronoi':
                nn_inputs = {'tol': 0.1, 'cutoff': radius_max}

        voronoi = CrystalNN(**nn_inputs) if nn_finder == 'crystal' else VoronoiNN(**nn_inputs)  # pylint: disable=unexpected-keyword-arg

        sites = self.hubbard_structure.sites
        pymat = self.hubbard_structure.get_pymatgen_structure()
        pairs = self.get_interacting_pairs()
        max_num_of_neighs = 0

        for i, site in enumerate(sites):  # pylint: disable=too-many-nested-blocks
            if site.kind_name in pairs:
                neigh_species = voronoi.get_cn_dict(pymat, i)  # e.g. {'O': 4, 'S': 2, ...}
                max_num_of_neighs = max(max_num_of_neighs, np.sum(list(neigh_species.values())))

        return max_num_of_neighs


def initialize_hubbard_parameters(
    structure: StructureData,
    pairs: Dict[str, Tuple[str, float, float, Dict[str, str]]],
    nn_finder: str = 'crystal',
    nn_inputs: Union[Dict, None] = None,
    fold: bool = True,
    standardize: bool = False,
    radius_max: float = 7.0,
    thr: float = 1e-5,
    use_kinds: bool = True,
    **_,
) -> HubbardStructureData:
    """Initialize the on-site and intersite parameters using nearest neighbour finders.

    It peforms a nearest neighbour analysis (via pymatgen modules) to find the first inersite
    neighbours for all the onsite atoms. Only the atoms in the `pair`, and that are nearest
    neighbours from the analysis, are initialized.

    :param structure: a StructureData instance
    :param pairs: dictionary of the kind
        {onsite name: [onsite manifold, onsite value, intersites value, {neighbour name: neighbour manifold}], ...}
        For example: {'Fe': ['3d', 5.0, 1.0, {'O':'2p', 'O1':'1s', 'Se':'4p'}]}
    :param nn_finder: string defining the nearest neighbour finder; options are:
        * `crystal`: use :class:`pymatgen.analysis.local_env.CrystalNN`
        * `voronoi`: use :class:`pymatgen.analysis.local_env.VoronoiNN`
    :param nn_inputs: inputs for the nearest neighbours finder; when None, standard inputs
        are used to find geometric first neighbours (recommended)
    :param fold: whether to fold in within the cell the atoms
    :param standardize: whether to standardize the atoms and the cell via spglib (symmetry analysis)
    :param radius_max: max radius where to look for neighbouring atoms, in Angstrom
    :param thr: threshold to refold the atoms with crystal coordinates close to 1.0
    :param use_kinds: whether to use kinds for initializing the parameters; when False, it
        initializes all the ``Kinds`` matching the given specie
    :return: HubbardStructureData with initialized Hubbard parameters
    """
    from aiida.tools.data import spglib_tuple_to_structure, structure_to_spglib_tuple
    from pymatgen.analysis.local_env import CrystalNN, VoronoiNN
    from spglib import standardize_cell

    if nn_finder not in ['crystal', 'voronoi']:
        raise ValueError('`nn_finder` must be either `cyrstal` or `voronoi`')

    if nn_inputs is None:
        if nn_finder == 'crystal':
            nn_inputs = {'distance_cutoffs': None, 'x_diff_weight': 0, 'porous_adjustment': False}
        if nn_finder == 'voronoi':
            nn_inputs = {'tol': 0.1, 'cutoff': radius_max}

    voronoi = CrystalNN(**nn_inputs) if nn_finder == 'crystal' else VoronoiNN(**nn_inputs)  # pylint: disable=unexpected-keyword-arg

    if not standardize and not fold:
        hubbard_structure = HubbardStructureData.from_structure(structure=structure)

    if fold:
        cell, kind_info, kinds = structure_to_spglib_tuple(structure)
        cell = list(cell)
        cell[1] = cell[1] % 1.0
        cell[1] = np.where(np.abs(cell[1] - 1.0) < thr, 0.0, cell[1])

        folded_struccture = spglib_tuple_to_structure(cell, kind_info, kinds)
        hubbard_structure = HubbardStructureData.from_structure(structure=folded_struccture)

    if standardize:
        cell, kind_info, kinds = structure_to_spglib_tuple(structure)
        new_cell = standardize_cell(cell)
        new_structure = spglib_tuple_to_structure(new_cell, kind_info, kinds)
        hubbard_structure = HubbardStructureData.from_structure(structure=new_structure)

    sites = hubbard_structure.sites
    name_to_specie = {kind.name: kind.symbol for kind in hubbard_structure.kinds}
    pymat = hubbard_structure.get_pymatgen_structure()

    for i, site in enumerate(sites):  # pylint: disable=too-many-nested-blocks
        site_name = site.kind_name if use_kinds else name_to_specie[site.kind_name]
        if site_name in pairs:
            neigh_species = voronoi.get_cn_dict(pymat, i)  # e.g. {'O': 4, 'S': 2, ...}
            onsite = pairs[site.kind_name]
            hubbard_structure.append_hubbard_parameter(i, onsite[0], i, onsite[0], onsite[1])

            for neigh_name in onsite[3]:
                try:  # we want an API that allows for specifying neighbours that are not present
                    specie = name_to_specie[neigh_name]

                    count = 0
                    if specie in neigh_species:
                        _, neigh_indices, images, distances = pymat.get_neighbor_list(sites=[pymat[i]], r=radius_max)

                        sort = np.argsort(distances)
                        neigh_indices = neigh_indices[sort]
                        images = images[sort]

                        for index, image in zip(neigh_indices, images):
                            if pymat[index].specie.name == specie:
                                hubbard_structure.append_hubbard_parameter(
                                    atom_index=i,
                                    atom_manifold=onsite[0],
                                    neighbour_index=int(index), # otherwise the validator complains
                                    neighbour_manifold=onsite[3][neigh_name],
                                    value=onsite[2],
                                    translation=np.array(image, dtype=np.int64).tolist(),
                                    hubbard_type='V',
                                )
                                count += 1

                            if count >= neigh_species[specie]:
                                break
                except KeyError:
                    pass

    return hubbard_structure


def get_supercell_atomic_index(index: int, num_sites: int, translation: List[Tuple[int, int, int]]) -> int:
    """Return the atomic index in 3x3x3 supercell.

    :param index: atomic index in unit cell
    :param num_sites: number of sites in structure
    :param translation: (3,) shape list of int referring to the translated atom in the 3x3x3 supercell

    :returns: atomic index in supercell standardized with the QuantumESPRESSO loop
    """
    return index + QE_TRANSLATIONS.index(translation) * num_sites


def get_index_and_translation(index: int, num_sites: int) -> Tuple[int, List[Tuple[int, int, int]]]:
    """Return the atomic index in unitcell and the associated translation from a 3x3x3 QuantumESPRESSO supercell index.

    :param index: atomic index
    :param num_sites: number of sites in structure
    :returns: tuple (index, (3,) shape list of ints)
    """
    from math import floor

    number = floor(index / num_sites)  # associated supercell number
    return (index - num_sites * number, QE_TRANSLATIONS[number])


def get_hubbard_indices(hubbard: Hubbard) -> List[int]:
    """Return the set list of Hubbard indices."""
    atom_indices = {parameters.atom_index for parameters in hubbard.parameters}
    neigh_indices = {parameters.neighbour_index for parameters in hubbard.parameters}
    atom_indices.update(neigh_indices)
    return list(atom_indices)


def is_intersite_hubbard(hubbard: Hubbard) -> bool:
    """Return whether `Hubbard` contains intersite interactions (+V)."""
    couples = [(param.atom_index != param.neighbour_index) for param in hubbard.parameters]
    return any(couples)


def max_number_of_neighbours(intersites_list: List[Tuple[int, int]]) -> int:
    """Return the maximum number of neighbours found.

    .. note:: it assumes only one onsite parameter is defined per atom index,
        that means `intra` Hubbard interactions are not defined, such as `V Fe 3d Fe 2p 1 1 1.0`

    :param intersites_list: list of lists of shape (atom index, neigbours index)
    """
    from collections import Counter

    first_indices = np.array(intersites_list, dtype='object')[:, 0]
    counts = Counter(first_indices)
    return max(counts.values()) - 1  # assuming there is always 1 and only 1 onsite parameter
