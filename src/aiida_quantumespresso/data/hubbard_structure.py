import numpy as np
from aiida import orm
from aiida.orm import StructureData

from copy import deepcopy


class HubbardStructureData(StructureData):
    """Structure data containing code independent info on Hubbard parameters."""
    
    def __init__(self, structure: StructureData, hubbard_parameters: list = None, **kwargs):
        """Constructor of HubbardStructureData.
        
        :param structure: StructureData
        :param hubbard_parameters: list of Hubbard parameters following the conventions
        """
        kwargs['cell'] = structure.cell
        
        super().__init__(**kwargs)
        
        # Inizializing the class attributes.
        for kind in structure.kinds:
            self.append_kind(kind)
        
        for site in structure.sites:
            self.append_site(site)
        
        self.hubbard_projectors = 'ortho-atomic'
        for hubbard_parameter in hubbard_parameters:
            self.append_hubbard_parameter(*hubbard_parameter)
    
    
    def atomic_distance(self, a, b, n=[0,0,0], m=[0,0,0]):
        """Compute distance between atom a and b in ase, translated by lattice vectors usign n and m.
        
        :params a, b: atomic indecis in unitcell
        :params n, m: (3,) arrays, traslation vectors of atom a and b (resp.) 
        """
        # positions in Cartesian coordinates
        ase = self.get_ase()
        
        R0a = ase.positions[a]
        R0b = ase.positions[b]
        
        cell = ase.cell
        
        Ra  = R0a + np.dot(n, cell)
        Rb  = R0b + np.dot(m, cell)
        
        return np.linalg.norm(Ra-Rb)
    
    def find_atom_images(self, a, b, dab, thr=1e-5):
        """
        This finds the integer vector m which traslates b to have distance dab from a.
        It can have no solutions, as well as many solutions.
        
        This is useful to get the indices in the QuantumESPRESSO logic.
        """
        traslations = self._get_standard_traslations()
        ns = []
        
        for traslation in traslations:
            d0ab = self.atomic_distance(a, b, [0,0,0], traslation)
            
            if np.abs(dab-d0ab) < thr:
                ns.append(traslation)
        
        return ns
    
    def find_atom_in_unitcell(self, position, thr=1e-5):
        """
        This finds the index of the atom within the unitcell from its image position.
        """
        ase = self.get_ase()
        
        cell = ase.cell
        inv_cell = np.linalg.inv(cell)
        positions = ase.positions
        
        for index, uc_position in enumerate(positions):
            diff = (uc_position-np.array(position))
            traslation = -np.dot(diff, inv_cell)
            traslation_int = np.rint(traslation)
            if np.all(np.isclose(traslation, traslation_int, thr)):
                break
        
        return index, traslation_int
    
    def get_images_distance_and_traslation(self, atom_index_i, atom_index_j):
        """Returns the distance and traslations between atom i and the 3x3x3 supercell images of atom j in ascending order in distance."""
        traslations = self._get_standard_traslations()
        distances = []
        
        for traslation in traslations:
            distance_ij = self.atomic_distance(atom_index_i, atom_index_j, [0,0,0], traslation)
            distances.append(distance_ij)
        
        sorting = np.argsort(distances)
        traslations = [traslations[i] for i in sorting]
        distances.sort()
        
        return distances, traslations
        
    @property
    def hubbard_projectors(self):
        """Get the Hubbard projectors."""
        return self.get_attribute('hubbard_projectors')
    
    @hubbard_projectors.setter
    def hubbard_projectors(self, value):
        self.set_attribute('hubbard_projectors', deepcopy(value))
        
    @property
    def hubbard_parameters(self):
        """Get the Hubbard parameters."""
        try:
            the_hps = self.get_attribute('hubbard_parameters')
        except AttributeError:
            the_hps = []
        return the_hps
    
    def append_hubbard_parameter(
        self, 
        atom_index_i, 
        orbital_i,
        atom_index_j,
        orbital_j,
        hubbard_value,
        hubbard_type="dudarev",
        distance_ij=None,
        traslation=None,
        thr=1e-5
    ):
        """Append a Hubbard parameter.
        
        ... inputs explanation ...
        If distance_ij is None, the smallest distance in the 3x3x3 supercell is taken.
        It raises error if more than one ideantical distance (within thr) is found.
        """
        if traslation is None:
            distances, _ = self.get_images_distance_and_traslation(atom_index_i, atom_index_j)
            
            if distance_ij is None:
                distance = distances[0]
            else:
                distance = distance_ij
            
            traslations = self.find_atom_images(atom_index_i, atom_index_j, distance, thr)
            
            if len(traslations) > 1:
                raise ValueError("more than one atom within threshold has been found with the same distance. Decrease threshold or provide traslation vector.")
            if len(traslations) == 0:
                raise ValueError("no atom image has been found within threshold. Check your distance or increase the threshold")
            
            traslation = traslations[0]
        
        hp = [
            atom_index_i,
            orbital_i, 
            atom_index_j,
            orbital_j, 
            hubbard_value,
            traslation,
            hubbard_type
        ]
        
        self._append_hubbard_parameter(hp)
    
    def pop_hubbard_parameter(self, index):
        """Pop a Hubbard parameter from the list."""
        hubbard_parameters = deepcopy(self.hubbard_parameters)
        hubbard_parameters.pop(index)
        
        self._set_hubbard_parameters(hubbard_parameters)
        
    def clean_hubbard_parameters(self):
        """Clean all the Hubbard parameters from the list."""
        self._set_hubbard_parameters([])
    
    def _set_hubbard_parameters(self, hubbard_parameters):
        """Set the full list of Hubbard parameters"""
        self.set_attribute("hubbard_parameters", hubbard_parameters)
    
    def _append_hubbard_parameter(self, hp_list):
        """Append to the Hubbard parameters a new hp."""
        hubbard_parameters = deepcopy(self.hubbard_parameters)
        hubbard_parameters.append(hp_list)
        
        self._set_hubbard_parameters(hubbard_parameters)
        
    def _get_standard_traslations(self):
        """Get first neighbours traslation vectors, i.e. 3x3x3 supercells. It is `standardized` to the QunatumESPRESSO loop."""
        from itertools import product
        return list(product((-1,0,1), repeat=3))
        
    def _get_quantum_espresso_hubbard_info(self):
        """Get 3x3x3 supercell atomic indecis, positions and types."""
        ase = self.get_ase()
        
        nat = ase.get_global_number_of_atoms()
        atom = nat

        uc_positions = ase.get_positions()
        uc_symbols= ase.get_chemical_symbols()
        cell = ase.cell

        sc_indecis = [i+1 for i in range(nat)] # QuantumESPRESSO starts from 1 
        sc_positions = [pos for pos in uc_positions]
        sc_symbols = [sym for sym in uc_symbols]

        traslations = self._get_standard_traslations()

        for traslation in traslations:
            if traslation!=(0,0,0):
                for na in range(nat):
                    atom += 1
                    sc_indecis.append(atom)
                    sc_positions.append(uc_positions[na]+np.dot(traslation, cell))
                    sc_symbols.append(uc_symbols[na])
        
        return sc_indecis, sc_positions, sc_symbols
        
    def get_quantum_espresso_hubbard_card(self):
        """Get QuantumESPRESSO `HUBBARD` input card for `pw.x`."""
        hubbard_parameters = deepcopy(self.hubbard_parameters)
        
        sites = self.sites
        natoms = len(sites)
        card = f'HUBBARD ({self.hubbard_projectors})\n'
        sc_indecis, _, _ = self._get_quantum_espresso_hubbard_info()
        
        traslations = self._get_standard_traslations()
        a = traslations.pop(13)
        traslations.insert(0, a)
        
        for hp in hubbard_parameters:
            atom_i = sites[hp[0]].kind_name
            atom_j = sites[hp[2]].kind_name
            orbital_i = hp[1]
            orbital_j = hp[3]
            value = hp [4]
            hubbard_type = hp[6]
            traslation = hp[5]
            
            for i, t in enumerate(traslations):
                if list(t)==list(traslation):
                    base_index = i
                    break
            
            index_i = hp[0]+1
            index_j = sc_indecis[base_index*natoms+hp[2]]
            
            if hubbard_type == "dudarev":
                pre = 'V'
            # --------------------------------------------------------------------------
            # Need to implement the logic for Liechtenstein as well. A bit more involved.
            # .
            # .
            # --------------------------------------------------------------------------
            
            line = f'{pre}\t{atom_i}-{orbital_i}\t{atom_j}-{orbital_j}\t{index_i}\t{index_j}\t{value}'
            line += '\n'
            card += line
        
        return card
    
    def get_quantum_espresso_hubbard_parameters(self):
        """Get QuantumESPRESSO `parameters.in` data for `pw.x` in `str`."""
        hubbard_parameters = deepcopy(self.hubbard_parameters)
        
        sites = self.sites
        natoms = len(sites)
        sc_indecis, _, _ = self._get_quantum_espresso_hubbard_info()
        card = ' # Atom 1  Atom 2  Hubbard V (eV)\n'
        
        traslations = self._get_standard_traslations()
        a = traslations.pop(13)
        traslations.insert(0, a)
        
        for hp in hubbard_parameters:
            value = hp [4]
            traslation = hp[5]
            hubbard_type = hp[6]
            
            for i, t in enumerate(traslations):
                if list(t)==list(traslation):
                    base_index = i
                    break
            
            index_i = hp[0]+1
            index_j = sc_indecis[base_index*natoms+hp[2]]
            
            if not hubbard_type == "dudarev":
                raise ValueError("`parameters.in` can be produced only with `dudarev`")
            
            line = f'\t{index_i}\t{index_j}\t{value}'
            line += '\n'
            card += line
        
        return card
    
    def find_first_neighbours_to_atom(self, atomic_index, neighbours_name, number_of_neighbours=6, use_kinds=True):
        """Find the first n-th atom neighbours with name to a specific atom (specified with the index).
        It gives back the neighbours indecis and traslation vetors."""
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
            raise ValueError("Neighbours not found in structure")
    
        distances = []
        traslations = []
        indecis = []
        
        for neighbour_index in neighbours_indecis:
            distances_n, traslations_n = self.get_images_distance_and_traslation(atomic_index, neighbour_index)
            indecis += [neighbour_index for _ in distances_n]
            distances += distances_n
            traslations += traslations_n
            
        sorting = np.argsort(distances)
        
        indecis = [indecis[i] for i in sorting]
        traslations = [traslations[i] for i in sorting]
        
        return indecis[:number_of_neighbours], traslations[:number_of_neighbours]
    
    def find_first_neighbours(self, atomic_name, neighbours_name, number_of_neighbours=6, use_kinds=True):
        """Find the first n-th atom neighbours with name to atoms with name.
        It gives back the atomic and neighbours indecis and the neighbours traslation vetors for each atom found."""
        sites = self.sites
        atom_indecis = []
        neighbours_indecis = []
        neighbours_traslations = []
        
        for i, site in enumerate(sites):
            if use_kinds:
                name = site.kind_name
            else:
                name = site.get_ase(self.kinds).symbol
            if name == atomic_name:
                atom_indecis.append(i)
        
        if not atom_indecis:
            raise ValueError("Atom not found in structure")
        
        for atomic_index in atom_indecis:
            indecis, traslations = self.find_first_neighbours_to_atom(atomic_index, neighbours_name, number_of_neighbours, use_kinds)
            neighbours_indecis.append(indecis)
            neighbours_traslations.append(traslations)
        
        return atom_indecis, neighbours_indecis, neighbours_traslations

    def find_first_neighbours_to_atom_within_radius(self, atomic_index, neighbours_name, radius, use_kinds=True):
        """Find the first n-th atom neighbours with name to a specific atom (specified with the index).
        It gives back the neighbours indecis and traslation vetors.
        
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
        
        if  not neighbours_indecis:
            raise ValueError("Neighbours not found in structure")
    
        distances = []
        traslations = []
        indecis = []
        
        for neighbour_index in neighbours_indecis:
            distances_n, traslations_n = self.get_images_distance_and_traslation(atomic_index, neighbour_index)
            indecis += [neighbour_index for _ in distances_n]
            distances += distances_n
            traslations += traslations_n
            
        sorting = []
        for i, distance in enumerate(distances):
            if distance < radius:
                sorting.append(i)
        
        indecis = [indecis[i] for i in sorting]
        traslations = [traslations[i] for i in sorting]
        
        return indecis, traslations
    
    def find_first_neighbours_within_radius(self, atom_name, neighbours_name, radius, use_kinds=True):
        """Find the first n-th atom neighbours with name to atoms with name.
        It gives back the atomic and neighbours indecis and the neighbours traslation vetors for each atom found.
        
        :param radius: radius in Angstrom
        """
        sites = self.sites
        atom_indecis = []
        neighbours_indecis = []
        neighbours_traslations = []
        
        for i, site in enumerate(sites):
            if use_kinds:
                name = site.kind_name
            else:
                name = site.get_ase(self.kinds).symbol
            if name == atom_name:
                atom_indecis.append(i)
        
        if not atom_indecis:
            raise ValueError("Atom not found in structure")
        
        for atomic_index in atom_indecis:
            indecis, traslations = self.find_first_neighbours_to_atom_within_radius(atomic_index, neighbours_name, radius, use_kinds)
            neighbours_indecis.append(indecis)
            neighbours_traslations.append(traslations)
        
        return atom_indecis, neighbours_indecis, neighbours_traslations
    
    def initialize_intersites_hubbard(
        self, 
        atom_name, 
        atom_orbital,
        neighbours_name,
        neighbours_orbital, 
        value,
        number_of_neighbours=None,
        radius=None,
        use_kinds=True
    ):
        """Initialize and append intersite Hubbard values between an atom and its neighbours."""
        if radius is not None:
            atom_indecis, neighbours_indecis, neighbours_traslations = self.find_first_neighbours_within_radius(atom_name, neighbours_name, radius=radius, use_kinds=use_kinds)
        elif number_of_neighbours is not None:
            atom_indecis, neighbours_indecis, neighbours_traslations = self.find_first_neighbours(atom_name, neighbours_name, number_of_neighbours=number_of_neighbours, use_kinds=use_kinds)
        else:
            raise ValueError("Either `number_of_neighbours` or `radius` must be specified.")
        
        for atom_index, neighbours_index, neighbours_traslation in zip(atom_indecis, neighbours_indecis, neighbours_traslations):
            for neighbour_index, neighbour_traslation in zip(neighbours_index, neighbours_traslation):
                self.append_hubbard_parameter(
                    atom_index_i=atom_index,
                    orbital_i=atom_orbital,
                    atom_index_j=neighbour_index,
                    orbital_j=neighbours_orbital,
                    hubbard_value=value,
                    traslation=neighbour_traslation,
                    hubbard_type="dudarev"
                )
                
    def hubbard_parameters_to_new_structure(self, new_structure: orm.StructureData, thr=1e-5):
        """It produces the Hubbard parameters for a new structure (e.g. a supercell) from the local parameters.
        
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
            
            uc_0_traslation = hubbard_parameter[5]

            uc_i_index = hubbard_parameter[0]
            uc_j_index = hubbard_parameter[2]
            
            uc_i_position = uc_positions[uc_i_index]
            uc_j_position = uc_positions[uc_j_index] 

            # print(np.linalg.norm(uc_j_position+np.dot(uc_0_traslation, uc_cell)-uc_i_position))
            
            sc_i_indecis = []
            sc_j_index = []
            sc_i_traslations = []

            for i, position in enumerate(sc_positions):
                traslation = np.dot(position-uc_i_position, uc_cell_inv)
                traslation_int = np.rint(traslation)
                if np.all(np.isclose(traslation, traslation_int, thr)):
                    sc_i_traslations.append(deepcopy(traslation_int))
                    sc_i_indecis.append(i)
                
                traslation = np.dot(position-uc_j_position, uc_cell_inv)
                traslation_int = np.rint(traslation)
                if np.all(np.isclose(traslation, traslation_int, thr)):
                    uc_j_traslation = np.array(traslation_int)
                    sc_j_index = i
            
            sc_j_position = sc_positions[sc_j_index]

            for sc_i_index, sc_i_traslation in zip(sc_i_indecis, sc_i_traslations):
                j_position = sc_j_position + np.dot(sc_i_traslation -uc_j_traslation +uc_0_traslation, uc_cell)
                new_j_index, new_traslation = sc_hps.find_atom_in_unitcell(j_position)
                
                # i_position = sc_positions[sc_i_index]
                # print(sc_hps.find_atom_in_unitcell(j_position))
                # print(np.linalg.norm(i_position-j_position))
                
                sc_hubbard_parameter = [
                    sc_i_index,
                    hubbard_parameter[1], 
                    new_j_index,
                    hubbard_parameter[3], 
                    hubbard_parameter[4],
                    new_traslation,
                    hubbard_parameter[6],
                ]
                
                sc_hubbard_parameters.append(sc_hubbard_parameter)
        
        return sc_hubbard_parameters
    
    def parse_quantum_espresso_hubbard_parameters(self, singlefile: orm.SinglefileData):
        """Parse a `parameters.in/out` file associated to the current structure and stores the Hubbard parameters."""
        hubbard_data = []
        
        with singlefile.open() as file:
            lines = file.readlines()
            for line in lines:
                if line.strip().split()[0] != '#':
                    hubbard_data.append([x for x in line.strip().split()])
        
        sites = self.sites
        natoms = len(sites)
        
        traslations = self._get_standard_traslations()
        a = traslations.pop(13)
        traslations.insert(0, a)
        
        found_i = False
        found_j = False
        
        for index_i, index_j, value in hubbard_data:
            for i, traslation in enumerate(traslations):
                if int(index_i) -(i+1)*natoms < 0 and not found_i:
                    index_0i = int(index_i) -i*natoms -1
                    found_i = True
                if int(index_j) -(i+1)*natoms < 0 and not found_j:
                    index_0j = int(index_j) -i*natoms -1
                    traslation_0 = traslation
                    found_j = True
            self.append_hubbard_parameter(
                atom_index_i=index_0i,
                orbital_i="",
                atom_index_j=index_0j,
                orbital_j="",
                traslation=traslation_0,
                hubbard_value=float(value),
                hubbard_type="dudarev",
            )
            found_j = False
            found_i = False
            
            
            