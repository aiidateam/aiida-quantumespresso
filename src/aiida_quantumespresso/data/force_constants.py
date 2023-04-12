# -*- coding: utf-8 -*-
"""Sub class of `Data` to handle interatomic force constants produced by the Quantum ESPRESSO q2r.x code."""
from aiida.orm import SinglefileData
import numpy
from qe_tools import CONSTANTS


class ForceConstantsData(SinglefileData):
    """Class to handle interatomic force constants from the Quantum ESPRESSO q2r.x code."""

    def set_file(self, file, filename=None, **kwargs):
        """Add a file to the node, parse it and set the attributes found.

        :param file: absolute path to the file or a filelike object
        :param filename: specify filename to use (defaults to name of provided file).
        """
        # pylint: disable=redefined-builtin
        super().set_file(file, filename, **kwargs)

        # Parse the force constants file
        dictionary, _, _ = parse_q2r_force_constants_file(self.get_content().splitlines(), also_force_constants=False)

        # Add all other attributes found in the parsed dictionary
        for key, value in dictionary.items():
            self.base.attributes.set(key, value)

    @property
    def number_of_species(self):
        """Return the number of atom species.

        :return: a scalar
        """
        return self.base.attributes.get('number_of_species')

    @property
    def number_of_atoms(self):
        """Return the number of atoms.

        :return: a scalar
        """
        return self.base.attributes.get('number_of_atoms')

    @property
    def cell(self):
        """Return the crystal unit cell where rows are the crystal vectors.

        :return: a 3x3 numpy.array
        """
        return numpy.array(self.base.attributes.get('cell'))

    @property
    def atom_list(self):
        """Return the list of atoms.

        :return: a list of length-5 tuple (element name, element mass amu_ry, 3 coordinates in cartesian Angstrom)
        """
        return self.base.attributes.get('atom_list')

    @property
    def has_done_electric_field(self):
        """Return flag to indicate if dielectric tensor and effective charges were computed.

        :return: a boolean
        """
        return self.base.attributes.get('has_done_electric_field')

    @property
    def dielectric_tensor(self):
        """Return the dielectric tensor matrix.

        :return: a 3x3 tuple
        """
        return self.base.attributes.get('dielectric_tensor')

    @property
    def effective_charges_eu(self):
        """Return the effective charges for each atom.

        :return: a list of number_of_atoms elements, each being a 3x3 tuple
        """
        return self.base.attributes.get('effective_charges_eu')

    @property
    def qpoints_mesh(self):
        """Return the number of q-points in each direction.

        :return: a length-3 tuple
        """
        return tuple(self.base.attributes.get('qpoints_mesh'))


def parse_q2r_force_constants_file(lines, also_force_constants=False):
    """Parse the real-space interatomic force constants file from QE-Q2R.

    :param also_force_constants: True to parse the force constants as well

    :return parsed_data: dictionary with the following keywords:

    - number_of_species: number of atom species ('ntyp' in QE)
    - number_of_atoms: number of atoms ('nat' in QE)
    - cell: unit cell
    - atom_list: list with, for each atom in the cell, a length-5
      tuple of the form (element_name, mass_in_amu_ry, then 3 coordinates in
      cartesian & Angstroms)
    - has_done_electric_field: True if dielectric constants & effective
      charges were computed
    - dielectric_tensor: dielectric constant (3x3 matrix)
    - effective_charges_eu: effective charges (ntyp x 3 x 3 matrix)
    - qpoints_mesh: length-3 tuple with number of qpoints in each dimension
      of the reciprocal lattice
    - force_constants: the real-space force constants: array with 7 indices, of the kind
        C(mi1, mi2, mi3, ji1, ji2, na1, na2) with
        * (mi1, mi2, mi3): the supercell dimensions
        * (ji1, ji2): axis of the displacement of the two atoms (from 1 to 3)
        * (na1, na2): atom numbers in the cell.
    - warnings: a list of warnings

    :return force_constants: the real-space force constants: array with 7 indices, of the kind
        C(mi1, mi2, mi3, ji1, ji2, na1, na2) where:
        * (mi1, mi2, mi3): the supercell dimensions
        * (ji1, ji2): axis of the displacement of the two atoms (from 1 to 3)
        * (na1, na2): atom numbers in the cell.
    """
    # pylint: disable=too-many-statements,too-many-branches,too-many-nested-blocks

    parsed_data = {}
    warnings = []

    try:
        # read first line
        current_line = 0
        first_line = lines[current_line].split()
        ntyp = int(first_line[0])
        nat = int(first_line[1])
        ibrav = int(first_line[2])
        celldm = [float(c) for c in first_line[3:]]
        if len(celldm) != 6:
            warnings.append('Wrong length for celldm')
        if ibrav != 0:
            warnings.append(f'ibrav ({ibrav}) is not 0; q-points path for phonon dispersion might be wrong')
        if any(item != 0 for item in celldm[1:]):
            warnings.append('celldm[1:] are not all zero; only celldm[0] will be used')

        parsed_data['number_of_species'] = ntyp
        parsed_data['number_of_atoms'] = nat
        current_line += 1

        # read cell data
        cell = tuple(
            tuple(float(c) * celldm[0] * CONSTANTS.bohr_to_ang
                  for c in l.split())
            for l in lines[current_line:current_line + 3]
        )
        parsed_data['cell'] = cell
        current_line += 3

        # read atom types and masses
        atom_type_list = []
        for ityp in range(ntyp):
            line = lines[current_line].split("'")
            if int(line[0]) == ityp + 1:
                atom_type_list.append(tuple((line[1].strip(), float(line[2]))))
            current_line += 1

        # read each atom coordinates
        atom_list = []
        for _ in range(nat):
            line = [float(c) for c in lines[current_line].split()]
            ityp = int(line[1])
            if 0 < ityp < ntyp + 1:
                line[0] = atom_type_list[ityp - 1][0]  # string with element name
                line[1] = atom_type_list[ityp - 1][1]  # element mass in amu_ry
                # Convert atomic positions (in cartesian) from alat to Angstrom:
                line[2:] = [pos * celldm[0] * CONSTANTS.bohr_to_ang for pos in line[2:]]
            atom_list.append(tuple(line))
            current_line += 1

        parsed_data['atom_list'] = atom_list

        # read lrigid (flag for dielectric constant and effective charges
        has_done_electric_field = lines[current_line].split()[0] == 'T'
        parsed_data['has_done_electric_field'] = has_done_electric_field
        current_line += 1

        if has_done_electric_field:
            # read dielectric tensor
            dielectric_tensor = tuple(tuple(float(c) for c in l.split()) for l in lines[current_line:current_line + 3])
            current_line += 3
            effective_charges_eu = []
            for _ in range(nat):
                current_line += 1
                effective_charges_eu.append(
                    tuple(tuple(float(c) for c in l.split()) for l in lines[current_line:current_line + 3])
                )
                current_line += 3

            parsed_data['dielectric_tensor'] = dielectric_tensor
            parsed_data['effective_charges_eu'] = effective_charges_eu

        # read q-points mesh
        qpoints_mesh = tuple(int(c) for c in lines[current_line].split())
        current_line += 1
        parsed_data['qpoints_mesh'] = qpoints_mesh

        force_constants = ()
        if also_force_constants:
            # read force_constants
            force_constants = numpy.zeros(qpoints_mesh + (3, 3, nat, nat), dtype=float)
            for ji1 in range(3):
                for ji2 in range(3):
                    for na1 in range(nat):
                        for na2 in range(nat):

                            indices = tuple(int(c) for c in lines[current_line].split())
                            current_line += 1
                            if (ji1 + 1, ji2 + 1, na1 + 1, na2 + 1) != indices:
                                raise ValueError('Wrong indices in force constants')

                            for mi3 in range(qpoints_mesh[2]):
                                for mi2 in range(qpoints_mesh[1]):
                                    for mi1 in range(qpoints_mesh[0]):

                                        line = lines[current_line].split()
                                        indices = tuple(int(c) for c in line[:3])

                                        if (mi1 + 1, mi2 + 1, mi3 + 1) != indices:
                                            raise ValueError('Wrong supercell indices in force constants')

                                        force_constants[mi1, mi2, mi3, ji1, ji2, na1, na2] = float(line[3])
                                        current_line += 1

    except (IndexError, ValueError) as exc:
        raise ValueError(str(exc) + '\nForce constants file could not be parsed (incorrect file format)') from exc

    return parsed_data, force_constants, warnings
