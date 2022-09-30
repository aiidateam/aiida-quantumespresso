# -*- coding: utf-8 -*-
from aiida.orm import Dict, TrajectoryData
import numpy
from packaging.version import Version
from qe_tools import CONSTANTS

from .base import Parser
from .parse_raw.cp import parse_cp_raw_output, parse_cp_traj_stanzas


class CpParser(Parser):
    """This class is the implementation of the Parser class for Cp."""

    def parse(self, **kwargs):
        """Receives in input a dictionary of retrieved nodes.

        Does all the logic here.
        """
        retrieved = self.retrieved

        # check what is inside the folder
        list_of_files = retrieved.base.repository.list_object_names()

        # options.metadata become attributes like this:
        stdout_filename = self.node.base.attributes.get('output_filename')
        # at least the stdout should exist
        if stdout_filename not in list_of_files:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_READ)

        # This should match 1 file
        xml_files = [xml_file for xml_file in self.node.process_class.xml_filenames if xml_file in list_of_files]
        if not xml_files:
            return self.exit(self.exit_codes.ERROR_MISSING_XML_FILE)
        elif len(xml_files) > 1:
            return self.exit(self.exit_codes.ERROR_OUTPUT_XML_MULTIPLE)

        # cp.x can produce, depending on the particular version of the code, a file called `print_counter.xml` or
        # `print_counter`, which is a plain text file with the number of the last timestep written in the trajectory
        # output. Note that if no trajectory is produced (for example because a single conjugate gradient step was
        # performed to calculate the ground state and the wavefunctions velocities) no printer_counter* file is written.

        print_counter_xml = True
        no_trajectory_output = False

        filename_counter_txt = self.node.process_class._FILE_PRINT_COUNTER_BASENAME
        filename_counter_xml = self.node.process_class._FILE_XML_PRINT_COUNTER_BASENAME

        # The following can happen and is not an error!
        if filename_counter_xml not in list_of_files and filename_counter_txt not in list_of_files:
            self.logger.error(
                f'We could not find the print counter file (`{filename_counter_txt}` or `{filename_counter_xml}`), '
                'assuming no trajectory output was produced'
            )
            no_trajectory_output = True

        if not no_trajectory_output:
            if filename_counter_txt in list_of_files:
                self.logger.info('print counter not in xml format')
                print_counter_xml = False
                filename_counter = filename_counter_txt
            else:  # xml format
                print_counter_xml = True
                self.logger.info('print counter in xml format')
                filename_counter = filename_counter_xml

        output_stdout = retrieved.base.repository.get_object_content(stdout_filename)
        output_xml = retrieved.base.repository.get_object_content(xml_files[0])
        output_xml_counter = None if no_trajectory_output else retrieved.base.repository.get_object_content(filename_counter)
        out_dict, _raw_successful = parse_cp_raw_output(
            output_stdout, output_xml, output_xml_counter, print_counter_xml
        )

        if not no_trajectory_output:
            # parse the trajectory. Units in Angstrom, picoseconds and eV.
            # append everthing in the temporary dictionary raw_trajectory
            raw_trajectory = {}
            evp_keys = [
                'electronic_kinetic_energy', 'cell_temperature', 'ionic_temperature', 'scf_total_energy', 'enthalpy',
                'enthalpy_plus_kinetic', 'energy_constant_motion', 'volume', 'pressure'
            ]

            # order of atom in the output trajectory changed somewhere after 6.5
            if Version(out_dict['creator_version']) > Version('6.5'):
                new_cp_ordering = True
            else:
                new_cp_ordering = False

            # Now prepare the reordering, as files in the xml are ordered
            if new_cp_ordering:
                reordering = None
            else:
                try:
                    # this works for old xml only
                    reordering = self._generate_sites_ordering(out_dict['species'], out_dict['atoms'])
                except KeyError:
                    # this works for newer versions
                    reordering = self._generate_sites_ordering(
                        out_dict['structure']['species'], out_dict['structure']['atoms']
                    )

            pos_filename = f'{self.node.process_class._PREFIX}.pos'
            if pos_filename not in list_of_files:
                out_dict['warnings'].append('Unable to open the POS file... skipping.')
                return self.exit_codes.ERROR_READING_POS_FILE
            number_of_atoms = out_dict.get(
                'number_of_atoms', out_dict['structure']['number_of_atoms'] if 'structure' in out_dict else None
            )
            trajectories = [
                ('positions', 'pos', CONSTANTS.bohr_to_ang, number_of_atoms),
                ('cells', 'cel', CONSTANTS.bohr_to_ang, 3),
                ('velocities', 'vel', CONSTANTS.bohr_to_ang / (CONSTANTS.timeau_to_sec * 10**12), number_of_atoms),
                ('forces', 'for', CONSTANTS.hartree_to_ev / CONSTANTS.bohr_to_ang, number_of_atoms),
                ('stresses', 'str', 1.0, 3) #stress in GPa
            ]

            for name, extension, scale, elements in trajectories:
                try:
                    with retrieved.base.repository.open(f'{self.node.process_class._PREFIX}.{extension}') as datafile:
                        data = [l.split() for l in datafile]
                        # POSITIONS stored in angstrom
                    traj_data = parse_cp_traj_stanzas(
                        num_elements=elements, splitlines=data, prepend_name=f'{name}_traj', rescale=scale
                    )
                    # here initialize the dictionary.
                    if extension == 'cel':
                        # NOTE: the trajectory output has the cell matrix transposed!!
                        raw_trajectory['cells'] = numpy.array(traj_data['cells_traj_data']).transpose((0, 2, 1))
                    elif extension == 'str':
                        raw_trajectory['stresses'] = numpy.array(traj_data['stresses_traj_data'])
                    else:
                        raw_trajectory[f'{name}_ordered'] = self._get_reordered_array(
                            traj_data[f'{name}_traj_data'], reordering
                        )
                    if extension == 'pos':
                        raw_trajectory['traj_times'] = numpy.array(traj_data[f'{name}_traj_times'])
                except IOError:
                    out_dict['warnings'].append(f'Unable to open the {extension.upper()} file... skipping.')

            # =============== EVP trajectory ============================
            try:
                with retrieved.base.repository.open(f'{self._node.process_class._PREFIX}.evp') as handle:
                    matrix = numpy.genfromtxt(handle)
                # there might be a different format if the matrix has one row only
                try:
                    matrix.shape[1]
                except IndexError:
                    matrix = numpy.array(numpy.matrix(matrix))

                if Version(out_dict['creator_version']) > Version('5.1'):
                    # Between version 5.1 and 5.1.1, someone decided to change
                    # the .evp output format, without any way to know that this
                    # happened... SVN commit 11158.
                    # I here use the version number to parse, plus some
                    # heuristics to check that I'm doing the right thing
                    #print "New version"
                    raw_trajectory['steps'] = numpy.array(matrix[:, 0], dtype=int)
                    raw_trajectory['times'] = matrix[:, 1]  # TPS, ps
                    raw_trajectory['electronic_kinetic_energy'] = matrix[:, 2] * CONSTANTS.hartree_to_ev  # EKINC, eV
                    raw_trajectory['cell_temperature'] = matrix[:, 3]  # TEMPH, K
                    raw_trajectory['ionic_temperature'] = matrix[:, 4]  # TEMPP, K
                    raw_trajectory['scf_total_energy'] = matrix[:, 5] * CONSTANTS.hartree_to_ev  # ETOT, eV
                    raw_trajectory['enthalpy'] = matrix[:, 6] * CONSTANTS.hartree_to_ev  # ENTHAL, eV
                    raw_trajectory['enthalpy_plus_kinetic'] = matrix[:, 7] * CONSTANTS.hartree_to_ev  # ECONS, eV
                    raw_trajectory['energy_constant_motion'] = matrix[:, 8] * CONSTANTS.hartree_to_ev  # ECONT, eV
                    raw_trajectory['volume'] = matrix[:, 9] * (CONSTANTS.bohr_to_ang**3)  # volume, angstrom^3
                    raw_trajectory['pressure'] = matrix[:, 10]  # out_press, GPa
                else:
                    #print "Old version"
                    raw_trajectory['steps'] = numpy.array(matrix[:, 0], dtype=int)
                    raw_trajectory['electronic_kinetic_energy'] = matrix[:, 1] * CONSTANTS.hartree_to_ev  # EKINC, eV
                    raw_trajectory['cell_temperature'] = matrix[:, 2]  # TEMPH, K
                    raw_trajectory['ionic_temperature'] = matrix[:, 3]  # TEMPP, K
                    raw_trajectory['scf_total_energy'] = matrix[:, 4] * CONSTANTS.hartree_to_ev  # ETOT, eV
                    raw_trajectory['enthalpy'] = matrix[:, 5] * CONSTANTS.hartree_to_ev  # ENTHAL, eV
                    raw_trajectory['enthalpy_plus_kinetic'] = matrix[:, 6] * CONSTANTS.hartree_to_ev  # ECONS, eV
                    raw_trajectory['energy_constant_motion'] = matrix[:, 7] * CONSTANTS.hartree_to_ev  # ECONT, eV
                    raw_trajectory['volume'] = matrix[:, 8] * (CONSTANTS.bohr_to_ang**3)  # volume, angstrom^3
                    raw_trajectory['pressure'] = matrix[:, 9]  # out_press, GPa
                    raw_trajectory['times'] = matrix[:, 10]  # TPS, ps

                # Huristics to understand if it's correct.
                # A better heuristics could also try to fix possible issues
                # (in new versions of QE, it's possible to recompile it with
                # the __OLD_FORMAT flag to get back the old version format...)
                # but I won't do it, as there may be also other columns swapped.
                # Better to stop and ask the user to check what's going on.

                #work around for 100ps format bug
                mask = numpy.array(raw_trajectory['traj_times']) >= 0
                len_bugged = len(numpy.array(raw_trajectory['times'])[mask == False])
                len_ok = len(numpy.array(raw_trajectory['times'])[mask])
                if len_ok > 0:
                    max_time_difference = abs(
                        numpy.array(raw_trajectory['times'])[mask] - numpy.array(raw_trajectory['traj_times'])[mask]
                    ).max()
                else:
                    max_time_difference = 0.0

                if max_time_difference > 1.e-4 or (
                    len_bugged > 0 and numpy.array(raw_trajectory['times'])[mask == False].min() < 100.0
                ):  # It is typically ~1.e-7 due to roundoff errors
                    # If there is a large discrepancy
                    # it means there is something very weird going on...
                    return self.exit_codes.ERROR_READING_TRAJECTORY_DATA

                # keep both times array (that usually are duplicated)
                # so that the user can check them by himselves
                if len_bugged > 0:
                    out_dict['warnings'].append(
                        '100ps format bug detected: ignoring trajectory\'s printed time from 100ps on'
                    )
            except IOError:
                out_dict['warnings'].append('Unable to open the EVP file... skipping.')

            # get the symbols from the input
            # TODO: I should have kinds in TrajectoryData
            input_structure = self.node.inputs.structure
            raw_trajectory['symbols'] = [str(i.kind_name) for i in input_structure.sites]

            traj = TrajectoryData()
            traj.set_trajectory(
                stepids=raw_trajectory['steps'],
                cells=raw_trajectory['cells'],
                symbols=raw_trajectory['symbols'],
                positions=raw_trajectory['positions_ordered'],
                times=raw_trajectory['times'],
                velocities=raw_trajectory['velocities_ordered'],
            )

            # eventually set the forces
            try:
                traj.set_array('forces', raw_trajectory['forces_ordered'])
            except KeyError:
                out_dict['warnings'].append('failed to set forces')

            # eventually set the stress
            if 'stresses' in raw_trajectory:
                traj.set_array('stresses',raw_trajectory['stresses'])

            for this_name in evp_keys:
                try:
                    traj.set_array(this_name, raw_trajectory[this_name])
                except KeyError:
                    # Some columns may have not been parsed, skip
                    pass

            self.out('output_trajectory', traj)

        # Remove big dictionaries that would be redundant
        # For atoms and cell, there is a small possibility that nothing is parsed but then probably nothing moved.
        for key in [
            'atoms', 'cell', 'ions_positions_stau', 'ions_positions_svel', 'ions_positions_taui', 'atoms_index_list',
            'atoms_if_pos_list', 'ions_positions_force', 'bands', 'structure'
        ]:
            out_dict.pop(key, None)

        # convert the dictionary into an AiiDA object
        output_params = Dict(out_dict)
        self.out('output_parameters', output_params)

    def get_linkname_trajectory(self):
        """Returns the name of the link to the output_structure (None if not present)"""
        return 'output_trajectory'

    def _generate_sites_ordering(self, raw_species, raw_atoms):
        """take the positions of xml and from file.pos of the LAST step and compare them."""
        # Examples in the comments are for species [Ba, O, Ti]
        # and atoms [Ba, Ti, O, O, O]

        # Dictionary to associate the species name to the idx
        # Example: {'Ba': 1, 'O': 2, 'Ti': 3}
        species_dict = {name: idx for idx, name in zip(raw_species['index'], raw_species['type'])}
        # List of the indices of the specie associated to each atom,
        # in the order specified in input
        # Example: (1,3,2,2,2)
        atoms_species_idx = [species_dict[a[0]] for a in raw_atoms]
        # I also attach the current position; important to convert to a list
        # Otherwise the iterator can be looped on only once!
        # Example: ((0,1),(1,3),(2,2),(3,2),(4,2))
        ref_atom_list = list(enumerate(atoms_species_idx))
        new_order_tmp = []
        # I reorder the atoms, first by specie, then in their order
        # This is the order used in output by CP!!
        # Example: ((0,1),(2,2),(3,2),(4,2),(1,3))
        for specie_idx in sorted(raw_species['index']):
            for elem in ref_atom_list:
                if elem[1] == specie_idx:
                    new_order_tmp.append(elem)
        # This is the new order that is printed in CP:
        # e.g. reordering[2] is the index of the atom, in the input
        # list of atoms, that is printed in position 2 (0-based, so the
        # third atom) in the CP output files.
        # Example: [0,2,3,4,1]
        reordering = [_[0] for _ in new_order_tmp]
        # I now need the inverse reordering, to put back in place
        # from the output ordering to the input one!
        # Example: [0,4,1,2,3]
        # Because in the final list (Ba, O, O, O, Ti)
        # the first atom Ba in the input is atom 0 in the CP output (the first),
        # the second atom Ti in the input is atom 4 (the fifth) in the CP output,
        # and so on
        sorted_indexed_reordering = sorted([(_[1], _[0]) for _ in enumerate(reordering)])
        reordering_inverse = [_[1] for _ in sorted_indexed_reordering]
        return reordering_inverse

    def _get_reordered_list(self, origlist, reordering):
        """Given a list to reorder, a list of integer positions with the new order, return the reordered list."""
        return [origlist[e] for e in reordering]

    def _get_reordered_array(self, _input, reordering):
        if reordering is not None:
            return numpy.array([self._get_reordered_list(i, reordering) for i in _input])
        else:
            return numpy.array(_input)
