# -*- coding: utf-8 -*-
# pylint: disable=invalid-name,too-many-lines
"""Data plugin that represents a crystal structure with Hubbard parameters."""
from copy import deepcopy
import json
from typing import Union

from aiida.orm import SinglefileData, StructureData
import numpy as np


class HubbardStructureData(StructureData):
    """Structure data containing code independent info on Hubbard parameters.

    Todo
    ----
    * maybe put parts in one or more `Mixin`?
    * general clean up
    * finish typings and docstrings
    * add a `is_ordered_atoms` method?
    * add (py)tests
    * refine name of variables
    * decide upon final structure of hubbard parameters
    * understand how much to be restrictive in the input checks

    Notes
    -----
    * the class is initialised from a StructureData object - more convenient to let this free?
    * `hubbard_parameters` could be very long: shall we dump them into e.g. yaml files or similar?
    * the sorting algorithm is not the fastest, but usually not needed, so is it ok?

    """

    # pylint: disable=too-many-public-methods

    _hubbard_filename = 'hubbard_parameters.json'

    def __init__(
        self,
        structure: StructureData,
        hubbard_parameters: list = None,
        hubbard_projectors: str = 'ortho-atomic',
        **kwargs
    ):
        """Construct an instance.

        Examples for ``hubbard_parameters``:
        * Structure containing 1 Fe atom ==> [[0,'3d',0,'3d',5.0,[0,0,0],'dudarev']]
        * Structure containing 1 Fe and 1 O atoms:
            - [[0,'3d',0,'3d',5.0,[0,0,0],'dudarev'], [0,'3d',1,'2p',1.0,[0,0,0],'dudarev']]

        Note that many methods are provided to make this definition much easier using atom names/kind names.

        :param structure: a StructureData object
        :param hubbard_parameters: list of Hubbard parameters following the convention:
            * the list must contain onsite/intersite parameters as list of list: e.g. [[parameters 1], ...]
            * each onsite/intersite parameters must have the information in this precise order:
                1. Hubbard atom index (in the structure sites) for the first atom
                2. Hubbard atom orbital/manifold for the first atom
                3. Hubbard atom index (in the structure sites) for the second atom (same for onsite/traditional Hubbard)
                4. Hubbard atom orbital/manifold for the second atom
                5. Value of the Hubbard parameter in eV
                6. Translation vector (useful for intersite to specify neighbouring atoms)
                   outside the unitcell; for traditional Hubbard, this is [0,0,0])
                7. Hubbard formulation (e.g. `dudarev`) and type (e.g. `dudarev-j`)

        :param hubbard_projectors: Hubbard projector type that are used for the occupations
        """
        kwargs['cell'] = structure.cell

        super().__init__(**kwargs)

        # Inizializing the class attributes.
        for kind in structure.kinds:
            self.append_kind(kind)

        for site in structure.sites:
            self.append_site(site)

        # Default we select `ortho-atomic`
        self.hubbard_projectors = hubbard_projectors

        if hubbard_parameters is not None:
            self.hubbard_parameters = hubbard_parameters

    def atomic_distance(self, a: int, b: int, n: list = None, m: list = None):
        """Compute distance between atom a and b in ase, translated by lattice vectors usign n and m.

        :params a, b: atomic indecis in unitcell
        :params n, m: (3,) arrays, translation vectors of atom a and b (resp.)
        """
        n = n or [0, 0, 0]
        m = m or [0, 0, 0]
        # positions in Cartesian coordinates
        ase = self.get_ase()

        R0a = ase.positions[a]
        R0b = ase.positions[b]

        cell = ase.cell

        Ra = R0a + np.dot(n, cell)
        Rb = R0b + np.dot(m, cell)

        return np.linalg.norm(Ra - Rb)

    def find_atom_images(self, a: int, b: int, dab: float, thr: float = 1e-5) -> list:
        """Find the integer vector m which traslates b to have distance dab from a.

        It can have no solutions, as well as many solutions.
        This is useful to get the indices in the QuantumESPRESSO logic.
        """
        translations = self._get_standard_translations()
        ns = []

        for translation in translations:
            d0ab = self.atomic_distance(a, b, [0, 0, 0], translation)

            if np.abs(dab - d0ab) < thr:
                ns.append(list(translation))

        return ns

    def find_atom_in_unitcell(self, position: list, thr: float = 1e-5):
        """Find the index of the atom within the unitcell from its image position.

        :param position: (3,) array list describing the position of the atom in Cartesian coordinates.
        """
        ase = self.get_ase()

        cell = ase.cell
        inv_cell = np.linalg.inv(cell)
        positions = ase.positions

        if not positions.any():
            raise ValueError('the structure has no positions')

        index = 0

        for index, uc_position in enumerate(positions):
            diff = (uc_position - np.array(position))
            translation = -np.dot(diff, inv_cell)
            translation_int = np.rint(translation)
            if np.all(np.isclose(translation, translation_int, thr)):
                break
        else:
            raise ValueError('failed to find the atom in the unit cell.')

        return index, translation_int

    def get_images_distance_and_translation(self, atom_index_i: int, atom_index_j: int):
        """Return the distance and translations between atom i and the 3x3x3 supercell images of atom j.

        Results are ordered by distance in ascending order.

        :param atom_index_*: number indicating the atom * in the unitcell; it ranges between 0 and (number of atoms -1)
        """
        translations = self._get_standard_translations()
        distances = []

        for translation in translations:
            distance_ij = self.atomic_distance(atom_index_i, atom_index_j, [0, 0, 0], translation)
            distances.append(distance_ij)

        sorting = np.argsort(distances)
        translations = [translations[i] for i in sorting]
        distances.sort()

        return distances, translations

    @property
    def hubbard_projectors(self):
        """Get the Hubbard projectors.

        :return: string
        """
        return self.get_attribute('hubbard_projectors')

    @hubbard_projectors.setter
    def hubbard_projectors(self, value: str):
        """Set the Hubbard projectors."""
        if isinstance(value, str):
            self.base.attributes.set('hubbard_projectors', deepcopy(value))
        else:
            raise ValueError(f'type <{type(value)}> is not correct; only class <{type(str)}> is accepted')

    @property
    def hubbard_parameters(self):
        """Get the Hubbard parameters.

        :return: list of list with 7 values describing the Hubbard parameters.
        """
        if self._hubbard_filename in self.base.repository.list_object_names():
            with self.base.repository.open(self._hubbard_filename, mode='rb') as handle:
                return json.load(handle)

        return []

    @hubbard_parameters.setter
    def hubbard_parameters(self, hp_parameters: list):
        """Set the full Hubbard parameters from a list (overrides completely the current ones)."""
        self._set_hubbard_parameters(hp_parameters)

    def _set_hubbard_parameters(self, hubbard_parameters: list):
        """Check and set the full list of Hubbard parameters."""
        num_atoms = len(self.sites)
        if isinstance(hubbard_parameters, list):
            for parameters in hubbard_parameters:
                if not len(parameters) == 7:
                    raise ValueError('each line of the parameters must contain 7 values')
                for i in (0, 2):
                    if not isinstance(parameters[i], int):
                        raise ValueError(f'the index {i} element must be integer')
                    if parameters[i] >= num_atoms or parameters[i] < 0:
                        raise ValueError(f'the index {i} element must be integer')
                for i in (1, 3):
                    if not isinstance(parameters[i], str):
                        raise ValueError(f'the index {i} element must be a string')
                    if not len(parameters[i]) in (2, 5):
                        raise ValueError(f'the index {i} element must contain two or 5 characters only')
                    if len(parameters[i]) == 2:
                        if not parameters[i][0] in [str(_ + 1) for _ in range(6)]:
                            raise ValueError(f'the quantum number {parameters[i][0]} is not correct')
                        if not parameters[i][1] in ['s', 'p', 'd', 'f', 'h']:
                            raise ValueError(f'the manifold number {parameters[i][1]} is not correct')
                    if len(parameters[i]) == 5:
                        if not parameters[i][2] == '-':
                            raise ValueError(f'the separator {parameters[i][0]} is not allowed. Only `-`.')
                        if not parameters[i][3] in [str(_ + 1) for _ in range(6)]:
                            raise ValueError(f'the quantum number {parameters[i][0]} is not correct')
                        if not parameters[i][4] in ['s', 'p', 'd', 'f', 'h']:
                            raise ValueError(f'the manifold number {parameters[i][1]} is not correct')
                if not isinstance(parameters[4], float):
                    raise ValueError(f'the index {i} element must be a float, representing the Hubbard value in eV')
                if not isinstance(parameters[5], (list, np.ndarray)):
                    raise ValueError(f'the index {i} element must be a list, representing a translation vector')
                if not len(parameters[5]) == 3:
                    raise ValueError(
                        f'the index {i} element must contain 3 integers, representing a translation vector'
                    )
                if not isinstance(parameters[6], str):
                    raise ValueError(f'the index {i} element must be a string, representing the Hubbard formulation')
        else:
            raise ValueError('Hubbard parameters must be given as a list.')

        serialized = json.dumps(hubbard_parameters)
        self.base.repository.put_object_from_bytes(serialized.encode('utf-8'), 'hubbard_parameters.json')

    def append_hubbard_parameter(
        self,
        atom_index_i: int,
        manifold_i: str,
        atom_index_j: int,
        manifold_j: str,
        hubbard_value: float,
        distance_ij: float = None,
        translation: list = None,
        hubbard_type: str = 'dudarev',
        thr: float = 1e-5
    ):
        """Append a Hubbard parameter.

        ... inputs explanation ...
        If distance_ij is None, the smallest distance in the 3x3x3 supercell is taken.
        It raises error if more than one ideantical distance (within thr) is found.
        """
        if translation is None:
            distances, _ = self.get_images_distance_and_translation(atom_index_i, atom_index_j)

            if distance_ij is None:
                distance = distances[0]
            else:
                distance = distance_ij

            translation = self.find_atom_images(atom_index_i, atom_index_j, distance, thr)

            if len(translation) > 1:
                message = 'more than one atom within threshold has been found with the same distance. \
                    The following equivalent translations have been found:\
                    \n {translation} \n \
                    Decrease threshold or provide on of the above translation vector.'

                raise ValueError(message)
            if len(translation) == 0:
                raise ValueError(
                    'no atom image has been found within threshold. Check your distance or increase the threshold'
                )

            translation = translation[0]

        hp = [atom_index_i, manifold_i, atom_index_j, manifold_j, hubbard_value, translation, hubbard_type]

        self._append_hubbard_parameter(hp)

    def pop_hubbard_parameter(self, index: int):
        """Pop a Hubbard parameter from the list."""
        hubbard_parameters = deepcopy(self.hubbard_parameters)
        hubbard_parameters.pop(index)

        self._set_hubbard_parameters(hubbard_parameters)

    def clean_hubbard_parameters(self):
        """Clean all the Hubbard parameters from the list."""
        self._set_hubbard_parameters([])

    def _append_hubbard_parameter(self, hp_list: list):
        """Append to the Hubbard parameters a new hp."""
        hubbard_parameters = deepcopy(self.hubbard_parameters)
        hubbard_parameters.append(hp_list)

        self._set_hubbard_parameters(hubbard_parameters)

    @staticmethod
    def _get_standard_translations() -> Union[list, tuple]:
        """Get first neighbours translation vectors, i.e. 3x3x3 supercells.

        It is `standardized` to the QunatumESPRESSO loop.
        """
        from itertools import product
        return list(product((-1, 0, 1), repeat=3))

    def _get_quantum_espresso_hubbard_info(self) -> tuple:
        """Get 3x3x3 supercell atomic indecis, positions and types."""
        ase = self.get_ase()

        nat = ase.get_global_number_of_atoms()
        atom = nat

        # uc = `unitcell`
        uc_positions = ase.get_positions()
        uc_symbols = ase.get_chemical_symbols()
        cell = ase.cell

        # sc = `supercell`
        sc_indecis = [i + 1 for i in range(nat)]  # QuantumESPRESSO starts from 1
        sc_positions = list(uc_positions)
        sc_symbols = list(uc_symbols)

        translations = self._get_standard_translations()

        for translation in translations:
            if translation != (0, 0, 0):
                for na in range(nat):
                    atom += 1
                    sc_indecis.append(atom)
                    sc_positions.append(uc_positions[na] + np.dot(translation, cell))
                    sc_symbols.append(uc_symbols[na])

        return sc_indecis, sc_positions, sc_symbols

    def get_quantum_espresso_hubbard_card(self) -> str:
        """Get QuantumESPRESSO `HUBBARD` input card for `pw.x` for versions > v.7.1."""
        hubbard_parameters = self.hubbard_parameters
        sites = self.sites
        natoms = len(sites)

        lines = [f'HUBBARD\t{self.hubbard_projectors}\n']
        sc_indecis, _, _ = self._get_quantum_espresso_hubbard_info()

        translations = self._get_standard_translations()
        a = translations.pop(13)
        translations.insert(0, a)

        for hp in hubbard_parameters:
            atom_i = sites[hp[0]].kind_name
            atom_j = sites[hp[2]].kind_name
            manifold_i = hp[1]
            manifold_j = hp[3]
            value = hp[4]
            hubbard_type = hp[6]
            translation = hp[5]

            for i, t in enumerate(translations):
                if list(t) == list(translation):
                    base_index = i
                    break

            index_i = hp[0] + 1
            index_j = sc_indecis[base_index * natoms + hp[2]]

            if hubbard_type.lower().startswith('dudarev'):
                if hubbard_type[-1].lower() == 'j':
                    pre = 'J'
                elif atom_i == atom_j and manifold_i == manifold_j:
                    pre = 'U'
                else:
                    pre = 'V'

                if pre in ['U', 'J']:
                    line = f'{pre}\t{atom_i}-{manifold_i}\t{value}'
                else:
                    line = f'{pre}\t{atom_i}-{manifold_i}\t{atom_j}-{manifold_j}\t{index_i}\t{index_j}\t{value}'
            elif hubbard_type.lower().startswith('liechtenstein-'):
                pre = hubbard_type[-1].upper()
                line = f'{pre}\t{atom_i}-{manifold_i} \t{value}'
            else:
                raise ValueError(f'Hubbard type {hubbard_type} is not recognised.')

            line += '\n'
            if line not in lines:
                lines.append(line)

        return ' '.join(lines)

    def get_quantum_espresso_hubbard_parameters(self) -> str:
        """Get QuantumESPRESSO `parameters.in` data for `pw.x` in `str`."""
        hubbard_parameters = deepcopy(self.hubbard_parameters)

        sites = self.sites
        natoms = len(sites)
        sc_indecis, _, _ = self._get_quantum_espresso_hubbard_info()
        card = ' # Atom 1  Atom 2  Hubbard V (eV)\n'

        translations = self._get_standard_translations()
        a = translations.pop(13)
        translations.insert(0, a)

        for hp in hubbard_parameters:
            value = hp[4]
            translation = hp[5]
            hubbard_type = hp[6]

            for i, t in enumerate(translations):
                if list(t) == list(translation):
                    base_index = i
                    break

            index_i = hp[0] + 1
            index_j = sc_indecis[base_index * natoms + hp[2]]

            if not hubbard_type == 'dudarev':
                raise ValueError('`parameters.in` can be produced only with `dudarev`')

            line = f'\t{index_i}\t{index_j}\t{value}'
            line += '\n'
            card += line

        return card

    def _get_quantum_espresso_ordered_params(self) -> list:
        """Return rearranged Hubbard parameters for multiple Hubbard parameters on atoms and couple of atoms."""
        # pylint: disable=too-many-statements
        hubbard_parameters = deepcopy(self.hubbard_parameters)
        hp = np.array(hubbard_parameters, dtype='object')

        # We first sort the first column
        hp_argsort = np.argsort(np.array(hp[:, 0], dtype='int'))
        hubbard_parameters = [hubbard_parameters[i] for i in hp_argsort]

        # Second, we sort the third column maintaining the first column
        at_index = deepcopy(hubbard_parameters[0][0])  # leading running Hubbard atomic index
        sub_hp = []  # sub array to reorder
        first_index = 0  # this index determines the sub array to reorder
        last_index = len(hubbard_parameters) - 1

        for i, params in enumerate(hubbard_parameters):
            if params[0] == at_index and i != last_index:
                sub_hp.append(params)
            else:
                if i == last_index:
                    sub_hp.append(params)
                if len(sub_hp) > 1:
                    sub_hp_argsort = np.argsort(np.array(sub_hp, dtype='object')[:, 2])
                    hps = deepcopy(hubbard_parameters)
                    for k, l in enumerate(sub_hp_argsort):
                        hubbard_parameters[k + first_index] = hps[l + first_index]
                first_index = i
                at_index = deepcopy(hubbard_parameters[i][0])
                sub_hp = [params]

        # Third, we sort the column of the vectors
        at_index_1 = deepcopy(hubbard_parameters[0][0])  # leading running Hubbard atomic index
        at_index_2 = deepcopy(hubbard_parameters[0][2])  # leading running Hubbard atomic index
        sub_hp = []  # sub array to reorder
        first_index = 0  # this index determines the sub array to reorder

        for i, params in enumerate(hubbard_parameters):
            if params[0] == at_index_1 and params[2] == at_index_2 and i != last_index:
                try:
                    sub_hp = sub_hp.tolist()
                except AttributeError:
                    pass
                sub_hp.append(params)
            else:
                if i == last_index:
                    sub_hp.append(params)
                if len(sub_hp) > 1:
                    sub_hp = [shp[5] for shp in sub_hp]
                    sub_hp = np.array(sub_hp, dtype='int')
                    sub_hp_argsort = np.argsort(
                        sub_hp.view('i8,i8,i8'), order=['f0', 'f1', 'f2'], axis=0, kind='stable'
                    ).view('int').flatten()
                    hps = deepcopy(hubbard_parameters)
                    for k, l in enumerate(sub_hp_argsort):
                        hubbard_parameters[k + first_index] = hps[l + first_index]
                first_index = i
                at_index_1 = deepcopy(hubbard_parameters[i][0])
                at_index_2 = deepcopy(hubbard_parameters[i][2])
                sub_hp = [params]

        # Finally, we sort with the QE logic of `background` and `standard` orbitals
        # NOTE: this will work only if all the 4 combinations are given
        at_index_1 = deepcopy(hubbard_parameters[0][0])  # leading running Hubbard atomic index
        at_index_2 = deepcopy(hubbard_parameters[0][2])  # leading running Hubbard atomic index
        at_translation = deepcopy(hubbard_parameters[0][5])  # leading running Hubbard atomic index
        sub_hp = []  # sub array to reorder
        first_index = 0  # this index determines the sub array to reorder

        for i, params in enumerate(hubbard_parameters):
            if params[0] == at_index_1 and params[2] == at_index_2 and params[5] == at_translation and i != last_index:
                try:
                    sub_hp = sub_hp.tolist()
                except AttributeError:
                    pass
                sub_hp.append(params)
            else:
                if i == last_index:
                    sub_hp.append(params)
                if len(sub_hp) == 4:
                    sub_hp_copy = np.array(sub_hp, dtype='object')
                    # Ordering first column of orbitals with QE logic
                    sub_hp_num = np.array(get_orbital_number_array(sub_hp_copy[:, 1]), dtype='int')
                    sub_hp_argsort = np.argsort(sub_hp_num)[::-1]
                    sub_hp = np.array([sub_hp[j] for j in sub_hp_argsort], dtype='object')
                    # Ordering second column of orbitals with QE logic
                    sub_hp_num = np.array(get_orbital_number_array(sub_hp[:, 3]), dtype='int')
                    sub_hp_argsort_2 = order_qe_manifolds(sub_hp_num)
                    # Combining the ordering
                    sub_hp_argsort = [sub_hp_argsort[j] for j in sub_hp_argsort_2]
                    hps = deepcopy(hubbard_parameters)
                    for k, l in enumerate(sub_hp_argsort):
                        hubbard_parameters[k + first_index] = hps[l + first_index]
                first_index = i
                at_index_1 = hubbard_parameters[i][0]
                at_index_2 = hubbard_parameters[i][2]
                at_translation = hubbard_parameters[i][5]
                sub_hp = [params]

        return hubbard_parameters

    def find_first_neighbours_to_atom(
        self, atomic_index: int, neighbours_name: str, number_of_neighbours: int = 6, use_kinds: bool = True
    ) -> tuple:
        """Find the first n-th atom neighbours with name to a specific atom (specified with the index).

        It gives back the neighbours indecis and translation vetors.
        """
        sites = self.sites
        neighbours_indecis = []

        for i, site in enumerate(sites):
            if use_kinds:
                name = site.kind_name
            else:
                name = site.get_ase(self.kinds).symbol
            if name == neighbours_name:
                neighbours_indecis.append(i)

        if not neighbours_indecis:
            raise ValueError('Neighbours not found in structure')

        distances = []
        translations = []
        indecis = []

        for neighbour_index in neighbours_indecis:
            distances_n, translations_n = self.get_images_distance_and_translation(atomic_index, neighbour_index)
            indecis += [neighbour_index for _ in distances_n]
            distances += distances_n
            translations += translations_n

        sorting = np.argsort(distances)

        indecis = [indecis[i] for i in sorting]
        translations = [translations[i] for i in sorting]

        return indecis[:number_of_neighbours], translations[:number_of_neighbours]

    def find_first_neighbours(self, atomic_name, neighbours_name, number_of_neighbours=6, use_kinds=True):
        """Find the first n-th atom neighbours with name to atoms with name.

        It gives back the atomic and neighbours indecis and the neighbours translation vetors for each atom found.
        """
        sites = self.sites
        atom_indecis = []
        neighbours_indecis = []
        neighbours_translations = []

        for i, site in enumerate(sites):
            if use_kinds:
                name = site.kind_name
            else:
                name = site.get_ase(self.kinds).symbol
            if name == atomic_name:
                atom_indecis.append(i)

        if not atom_indecis:
            raise ValueError('Atom not found in structure')

        for atomic_index in atom_indecis:
            indecis, translations = self.find_first_neighbours_to_atom(
                atomic_index, neighbours_name, number_of_neighbours, use_kinds
            )
            neighbours_indecis.append(indecis)
            neighbours_translations.append(translations)

        return atom_indecis, neighbours_indecis, neighbours_translations

    def find_first_neighbours_to_atom_within_radius(self, atomic_index, neighbours_name, radius, use_kinds=True):
        """Find the first n-th atom neighbours with name to a specific atom (specified with the index).

        It gives back the neighbours indecis and translation vetors.

        :param radius: radius in Angstrom
        """
        sites = self.sites
        neighbours_indecis = []

        for i, site in enumerate(sites):
            if use_kinds:
                name = site.kind_name
            else:
                name = site.get_ase(self.kinds).symbol
            if name == neighbours_name:
                neighbours_indecis.append(i)

        if not neighbours_indecis:
            raise ValueError('Neighbours not found in structure')

        distances = []
        translations = []
        indecis = []

        for neighbour_index in neighbours_indecis:
            distances_n, translations_n = self.get_images_distance_and_translation(atomic_index, neighbour_index)
            indecis += [neighbour_index for _ in distances_n]
            distances += distances_n
            translations += translations_n

        sorting = []
        for i, distance in enumerate(distances):
            if distance < radius:
                sorting.append(i)

        indecis = [indecis[i] for i in sorting]
        translations = [translations[i] for i in sorting]

        return indecis, translations

    def find_first_neighbours_within_radius(self, atom_name, neighbours_name, radius, use_kinds=True):
        """Find the first n-th atom neighbours with name to atoms with name.

        It gives back the atomic and neighbours indecis and the neighbours translation vetors for each atom found.

        :param radius: radius in Angstrom
        """
        sites = self.sites
        atom_indecis = []
        neighbours_indecis = []
        neighbours_translations = []

        for i, site in enumerate(sites):
            if use_kinds:
                name = site.kind_name
            else:
                name = site.get_ase(self.kinds).symbol
            if name == atom_name:
                atom_indecis.append(i)

        if not atom_indecis:
            raise ValueError('Atom not found in structure')

        for atomic_index in atom_indecis:
            indecis, translations = self.find_first_neighbours_to_atom_within_radius(
                atomic_index, neighbours_name, radius, use_kinds
            )
            neighbours_indecis.append(indecis)
            neighbours_translations.append(translations)

        return atom_indecis, neighbours_indecis, neighbours_translations

    def initialize_intersites_hubbard(
        self,
        atom_name,
        atom_manifold,
        neighbours_name,
        neighbours_manifold,
        value,
        number_of_neighbours=None,
        radius=None,
        use_kinds=True,
        hubbard_type='dudarev'
    ):
        """Initialize and append intersite Hubbard values between an atom and its neighbours."""
        if radius is not None:
            atom_indecis, neighbours_indecis, neighbours_translations = self.find_first_neighbours_within_radius(
                atom_name, neighbours_name, radius=radius, use_kinds=use_kinds
            )
        elif number_of_neighbours is not None:
            atom_indecis, neighbours_indecis, neighbours_translations = self.find_first_neighbours(
                atom_name, neighbours_name, number_of_neighbours=number_of_neighbours, use_kinds=use_kinds
            )
        else:
            raise ValueError('Either `number_of_neighbours` or `radius` must be specified.')

        for atom_index, neighbours_index, neighbours_translation in zip(
            atom_indecis, neighbours_indecis, neighbours_translations
        ):
            for neighbour_index, neighbour_translation in zip(neighbours_index, neighbours_translation):
                self.append_hubbard_parameter(
                    atom_index_i=atom_index,
                    manifold_i=atom_manifold,
                    atom_index_j=neighbour_index,
                    manifold_j=neighbours_manifold,
                    hubbard_value=value,
                    translation=list(neighbour_translation),
                    hubbard_type=hubbard_type
                )

    def initialize_onsites_hubbard(self, atom_name, atom_manifold, value, use_kinds=True, hubbard_type='dudarev'):
        """Initialize and append onsite Hubbard values of atoms with specific name."""
        self.initialize_intersites_hubbard(
            atom_name,
            atom_manifold,
            atom_name,
            atom_manifold,
            value,
            number_of_neighbours=1,
            use_kinds=use_kinds,
            hubbard_type=hubbard_type
        )

    def hubbard_parameters_to_new_structure(self, new_structure: StructureData, thr=1e-5):
        """Return the Hubbard parameters for a new structure (e.g. a supercell) from the local parameters.

        .. note:: the two structure need to be commensurate within thr.
        """
        # positions in Cartesian coordinates
        uc_ase = self.get_ase()
        sc_ase = new_structure.get_ase()

        uc_positions = uc_ase.positions
        sc_positions = sc_ase.positions

        uc_cell = uc_ase.cell
        uc_cell_inv = np.linalg.inv(uc_cell)

        uc_hubbard_parameters = self.hubbard_parameters
        sc_hubbard_parameters = []

        sc_hps = HubbardStructureData(structure=new_structure)

        for hubbard_parameter in uc_hubbard_parameters:

            uc_0_translation = hubbard_parameter[5]

            uc_i_index = hubbard_parameter[0]
            uc_j_index = hubbard_parameter[2]

            uc_i_position = uc_positions[uc_i_index]
            uc_j_position = uc_positions[uc_j_index]

            sc_i_indecis = []
            sc_j_index = []
            sc_i_translations = []

            for i, position in enumerate(sc_positions):
                translation = np.dot(position - uc_i_position, uc_cell_inv)
                translation_int = np.rint(translation)
                if np.all(np.isclose(translation, translation_int, thr)):
                    sc_i_translations.append(deepcopy(translation_int))
                    sc_i_indecis.append(i)

                translation = np.dot(position - uc_j_position, uc_cell_inv)
                translation_int = np.rint(translation)
                if np.all(np.isclose(translation, translation_int, thr)):
                    uc_j_translation = np.array(translation_int)
                    sc_j_index = i

            sc_j_position = sc_positions[sc_j_index]

            for sc_i_index, sc_i_translation in zip(sc_i_indecis, sc_i_translations):
                j_position = sc_j_position + np.dot(sc_i_translation - uc_j_translation + uc_0_translation, uc_cell)
                new_j_index, new_translation = sc_hps.find_atom_in_unitcell(j_position)

                sc_hubbard_parameter = [
                    sc_i_index,
                    hubbard_parameter[1],
                    new_j_index,
                    hubbard_parameter[3],
                    hubbard_parameter[4],
                    new_translation.tolist(),
                    hubbard_parameter[6],
                ]

                sc_hubbard_parameters.append(sc_hubbard_parameter)

        return sc_hubbard_parameters

    def parse_quantum_espresso_hubbard_parameters(self, singlefile: SinglefileData):
        """Parse a `parameters.in/out` file associated to the current structure and stores the Hubbard parameters."""
        hubbard_data = []

        with singlefile.open() as file:
            lines = file.readlines()
            for line in lines:
                if line.strip().split()[0] != '#':
                    hubbard_data.append(list(line.strip().split()))

        sites = self.sites
        natoms = len(sites)

        translations = self._get_standard_translations()
        a = translations.pop(13)
        translations.insert(0, a)

        found_i = False
        found_j = False

        for index_i, index_j, value in hubbard_data:
            for i, translation in enumerate(translations):
                if int(index_i) - (i + 1) * natoms <= 0 and not found_i:
                    index_0i = int(index_i) - i * natoms - 1
                    found_i = True
                if int(index_j) - (i + 1) * natoms <= 0 and not found_j:
                    index_0j = int(index_j) - i * natoms - 1
                    translation_0 = translation
                    found_j = True
            self.append_hubbard_parameter(
                atom_index_i=index_0i,
                manifold_i='',
                atom_index_j=index_0j,
                manifold_j='',
                translation=translation_0,
                hubbard_value=float(value),
                hubbard_type='dudarev',
            )
            found_j = False
            found_i = False

    def parse_quantum_espresso_hubbard_dat(self, filepath):
        """Parse the `HUBBARD.dat` of QuantumESPRESSO v.>7.1 file associated to the current structure.

        Also stores the parsed Hubbard parameters.
        """
        # pylint: disable=too-many-statements
        hubbard_data = []

        with open(filepath, encoding='utf-8') as file:
            lines = file.readlines()
            for line in lines:
                if line.strip().split()[0] != '#':
                    hubbard_data.append(list(line.strip().split()))

        sites = self.sites
        natoms = len(sites)

        translations = self._get_standard_translations()
        a = translations.pop(13)
        translations.insert(0, a)

        found_i = False
        found_j = False

        hubbard_projectors = hubbard_data.pop(0)[1]
        if hubbard_projectors.startswith('('):
            hubbard_projectors = hubbard_projectors[1:-1]

        self.hubbard_projectors = hubbard_projectors

        for data in hubbard_data:
            if data[0] == 'V':
                index_i = int(data[3])
                index_j = int(data[4])

                # Brute force extraction of manifolds strings
                manifold_i = ''
                manifold_j = ''
                mi = data[1].strip().split('-')
                mj = data[2].strip().split('-')
                mi.pop(0)
                mj.pop(0)
                for m in mi:
                    manifold_i += f'{m}'
                for m in mj:
                    manifold_j += f'{m}'
                if len(mi) > 1:
                    manifold_i = manifold_i[:-1]
                if len(mj) < 1:
                    manifold_j = manifold_j[:-1]

                for i, translation in enumerate(translations):
                    if index_i - (i + 1) * natoms <= 0 and not found_i:
                        index_0i = index_i - i * natoms - 1
                        found_i = True
                    if index_j - (i + 1) * natoms <= 0 and not found_j:
                        index_0j = index_j - i * natoms - 1
                        translation_0 = translation
                        found_j = True

                value = float(data[5])

                self.append_hubbard_parameter(
                    atom_index_i=index_0i,
                    manifold_i=manifold_i,
                    atom_index_j=index_0j,
                    manifold_j=manifold_j,
                    hubbard_value=value,
                    translation=list(translation_0),
                    hubbard_type='dudarev',
                )
                found_j = False
                found_i = False
            else:
                manifold = data[1].strip().split('-')  # with name in front
                name = manifold.pop(0)
                # Brute force extraction of manifolds strings
                manifold_complete = ''
                for m in manifold:
                    manifold_complete += f'{m}'
                if len(manifold) > 1:
                    manifold_complete = manifold_complete[:-1]
                value = float(data[2])

                self.initialize_onsites_hubbard(name, manifold, value)

    def reorder_atoms(self):
        """Reorder the atoms following the QuantumESPRESSO standards, meaning Hubbard atoms must come first.

        .. note:: this will not work if the data node has been already stored!
        """
        new_structure = deepcopy(self)
        new_structure.clear_kinds()

        hub_param = np.array(self.hubbard_parameters, dtype=object)

        sites = deepcopy(self.sites)
        indices = hub_param[:, 0].tolist() + hub_param[:, 2].tolist()
        indices = list(set(indices))
        hubbard_kinds = []

        for index in indices:
            hubbard_kinds.append(sites[index].kind_name)

        hubbard_kinds = list(set(hubbard_kinds))
        hubbard_kinds.sort(reverse=True)

        ordered_sites = []

        while hubbard_kinds:

            hubbard_kind = hubbard_kinds.pop()

            hubbard_sites = []
            remaining_sites = []

            hubbard_sites = [s for s in sites if s.kind_name == hubbard_kind]
            remaining_sites = [s for s in sites if not s.kind_name == hubbard_kind]

            ordered_sites.extend(hubbard_sites)
            sites = deepcopy(remaining_sites)

        # Extend the current site list with the remaining non-hubbard sites
        ordered_sites.extend(sites)

        for site in ordered_sites:

            if site.kind_name not in new_structure.get_kind_names():
                kind = self.get_kind(site.kind_name)
                new_structure.append_kind(kind)

            new_structure.append_site(site)

        new_hubbard_parameters = self.hubbard_parameters_to_new_structure(new_structure=new_structure)

        self.clear_kinds()

        for kind in new_structure.kinds:
            self.append_kind(kind)
        for site in new_structure.sites:
            self.append_site(site)

        self.clean_hubbard_parameters()
        self.hubbard_parameters = new_hubbard_parameters


def get_orbital_number_array(array: list) -> list:
    """Return an array with conventional quantum letters into an array with corresponding integers."""
    return [assign_orbital_number(orbital_letter[-1]) for orbital_letter in array]


def assign_orbital_number(orbital_letter: str) -> int:
    """Return the  angular quantum number from letter."""
    if orbital_letter == 's':
        return 1
    if orbital_letter == 'p':
        return 2
    if orbital_letter == 'd':
        return 3
    if orbital_letter == 'f':
        return 4
    if orbital_letter == 'h':
        return 5

    raise ValueError(f'{orbital_letter} is not an angular quantum letter.')


def order_qe_manifolds(array: list) -> list:
    """Return the argsort order of the manifolds of the second column of `background` and `standard` manifolds.

    .. note: we are expecting an array of dimension (4,). We simply compare the first half and the second half.
        The order should be [greater_1, lower_1, lower_2, greater_2], where 1 and 2 refers to the frist and second half
        respectively.
    """
    sorting = [0, 1, 2, 3]

    if array[0] < array[1]:
        sorting[0] = 1
        sorting[1] = 0
    if array[2] > array[3]:
        sorting[2] = 3
        sorting[3] = 2

    return sorting
