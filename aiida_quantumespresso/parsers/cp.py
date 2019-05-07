# -*- coding: utf-8 -*-
from __future__ import absolute_import

from distutils.version import LooseVersion

import numpy
from aiida.common import NotExistent
from aiida.orm import Dict, TrajectoryData
from aiida.parsers import Parser
from six.moves import zip

from aiida_quantumespresso.parsers.constants import (bohr_to_ang,
                                                     hartree_to_ev,
                                                     timeau_to_sec)
from aiida_quantumespresso.parsers.raw_parser_cp import (parse_cp_raw_output,
                                                         parse_cp_traj_stanzas)


class CpParser(Parser):
    """
    This class is the implementation of the Parser class for Cp.
    """

    def parse(self, **kwargs):
        """
        Receives in input a dictionary of retrieved nodes.
        Does all the logic here.
        """
        try:
            out_folder = self.retrieved
        except NotExistent:
            self.logger.error("No retrieved folder found")
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # check what is inside the folder
        list_of_files = out_folder._repository.list_object_names()

        # options.metadata become attributes like this:
        stdout_filename = self.node.get_attribute('output_filename')
        # at least the stdout should exist
        if stdout_filename not in list_of_files:
            self.logger.error("Standard output not found")
            return self.exit_codes.ERROR_READING_OUTPUT_FILE

        # This should match 1 file
        xml_files = [
            xml_file for xml_file in self.node.process_class.xml_filenames
            if xml_file in list_of_files
        ]
        if not xml_files:
            self.logger.error('no XML output files found, which is required for parsing')
            return self.exit_codes.ERROR_MISSING_XML_FILE
        elif len(xml_files) > 1:
            self.logger.error('more than one XML file retrieved, which should never happen')
            return self.exit_codes.ERROR_MULTIPLE_XML_FILES

        if self.node.process_class._FILE_XML_PRINT_COUNTER_BASENAME not in list_of_files:
            self.logger.error('We could not find the print counter file in the output')
            # TODO: Add an error for this counter
            return self.exit_codes.ERROR_MISSING_XML_FILE

        # Let's pass file handlers to this function
        out_dict, _raw_successful = parse_cp_raw_output(
            out_folder.open(stdout_filename),
            out_folder.open(xml_files[0]),
            out_folder.open(self.node.process_class._FILE_XML_PRINT_COUNTER_BASENAME)
        )

        # parse the trajectory. Units in Angstrom, picoseconds and eV.
        # append everthing in the temporary dictionary raw_trajectory
        raw_trajectory = {}
        evp_keys = [
            'electronic_kinetic_energy', 'cell_temperature', 'ionic_temperature',
            'scf_total_energy', 'enthalpy', 'enthalpy_plus_kinetic',
            'energy_constant_motion', 'volume', 'pressure'
        ]

        # Now prepare the reordering, as filex in the xml are  ordered
        reordering = self._generate_sites_ordering(out_dict['species'],
                                                   out_dict['atoms'])

        pos_filename = '{}.{}'.format(self.node.process_class._PREFIX, 'pos')
        if pos_filename not in list_of_files:
            out_dict['warnings'].append("Unable to open the POS file... skipping.")
            return self.exit_codes.ERROR_READING_POS_FILE

        trajectories = [
            ('positions', 'pos', bohr_to_ang, out_dict['number_of_atoms']),
            ('cells', 'cel', bohr_to_ang, 3),
            ('velocities', 'vel', bohr_to_ang / timeau_to_sec * 10 ** 12, out_dict['number_of_atoms']),
        ]

        for name, extension, scale, elements in trajectories:
            try:
                with out_folder.open('{}.{}'.format(self.node.process_class._PREFIX, extension)) as datafile:
                    data = [l.split() for l in datafile]
                    # POSITIONS stored in angstrom
                traj_data = parse_cp_traj_stanzas(
                    num_elements=elements,
                    splitlines=data,
                    prepend_name='{}_traj'.format(name),
                    rescale=scale
                )
                # here initialize the dictionary. If the parsing of positions fails, though, I don't have anything
                # out of the CP dynamics. Therefore, the calculation status is set to FAILED.
                if extension != 'cel':
                    raw_trajectory['{}_ordered'.format(name)] = self._get_reordered_array(
                        traj_data['{}_traj_data'.format(name)],
                        reordering
                    )
                else:
                    raw_trajectory['cells'] = numpy.array(traj_data['cells_traj_data'])
                if extension == 'pos':
                    raw_trajectory['times'] = numpy.array(traj_data['{}_traj_times'.format(name)])
            except IOError:
                out_dict['warnings'].append("Unable to open the {} file... skipping.".format(extension.upper()))

        # =============== EVP trajectory ============================
        try:
            matrix = numpy.genfromtxt(out_folder.open('{}.evp'.format(self._node.process_class._PREFIX)))
            # there might be a different format if the matrix has one row only
            try:
                matrix.shape[1]
            except IndexError:
                matrix = numpy.array(numpy.matrix(matrix))

            if LooseVersion(out_dict['creator_version']) > LooseVersion("5.1"):
                # Between version 5.1 and 5.1.1, someone decided to change
                # the .evp output format, without any way to know that this
                # happened... SVN commit 11158.
                # I here use the version number to parse, plus some
                # heuristics to check that I'm doing the right thing
                #print "New version"
                raw_trajectory['steps'] = numpy.array(matrix[:,0],dtype=int)
                raw_trajectory['evp_times']                 = matrix[:,1]                    # TPS, ps
                raw_trajectory['electronic_kinetic_energy'] = matrix[:,2] * hartree_to_ev    # EKINC, eV
                raw_trajectory['cell_temperature']          = matrix[:,3]                    # TEMPH, K
                raw_trajectory['ionic_temperature']         = matrix[:,4]                    # TEMPP, K
                raw_trajectory['scf_total_energy']          = matrix[:,5] * hartree_to_ev    # ETOT, eV
                raw_trajectory['enthalpy']                  = matrix[:,6] * hartree_to_ev    # ENTHAL, eV
                raw_trajectory['enthalpy_plus_kinetic']     = matrix[:,7] * hartree_to_ev    # ECONS, eV
                raw_trajectory['energy_constant_motion']    = matrix[:,8] * hartree_to_ev    # ECONT, eV
                raw_trajectory['volume']                    = matrix[:,9] * (bohr_to_ang**3) # volume, angstrom^3
                raw_trajectory['pressure']                  = matrix[:,10]                    # out_press, GPa
            else:
                #print "Old version"
                raw_trajectory['steps'] = numpy.array(matrix[:,0],dtype=int)
                raw_trajectory['electronic_kinetic_energy'] = matrix[:,1] * hartree_to_ev    # EKINC, eV
                raw_trajectory['cell_temperature']          = matrix[:,2]                    # TEMPH, K
                raw_trajectory['ionic_temperature']         = matrix[:,3]                    # TEMPP, K
                raw_trajectory['scf_total_energy']          = matrix[:,4] * hartree_to_ev    # ETOT, eV
                raw_trajectory['enthalpy']                  = matrix[:,5] * hartree_to_ev    # ENTHAL, eV
                raw_trajectory['enthalpy_plus_kinetic']     = matrix[:,6] * hartree_to_ev    # ECONS, eV
                raw_trajectory['energy_constant_motion']    = matrix[:,7] * hartree_to_ev    # ECONT, eV
                raw_trajectory['volume']                    = matrix[:,8] * (bohr_to_ang**3) # volume, angstrom^3
                raw_trajectory['pressure']                  = matrix[:,9]                    # out_press, GPa
                raw_trajectory['evp_times']                  = matrix[:,10]                    # TPS, ps

            # Huristics to understand if it's correct.
            # A better heuristics could also try to fix possible issues
            # (in new versions of QE, it's possible to recompile it with
            # the __OLD_FORMAT flag to get back the old version format...)
            # but I won't do it, as there may be also other columns swapped.
            # Better to stop and ask the user to check what's going on.
            max_time_difference = abs(
                numpy.array(raw_trajectory['times']) -
                numpy.array(raw_trajectory['evp_times'])).max()
            if max_time_difference > 1.e-4: # It is typically ~1.e-7 due to roundoff errors
                # If there is a large discrepancy
                # it means there is something very weird going on...
                return self.exit_codes.ERROR_READING_TRAJECTORY_DATA

            # Delete evp_times in any case, it's a duplicate of 'times'
            del raw_trajectory['evp_times']
        except IOError:
            out_dict['warnings'].append("Unable to open the EVP file... skipping.")

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

        for this_name in evp_keys:
            try:
                traj.set_array(this_name,raw_trajectory[this_name])
            except KeyError:
                # Some columns may have not been parsed, skip
                pass

        self.out('output_trajectory', traj)

        # Remove big dictionaries that would be redundant
        # For atoms and cell, there is a small possibility that nothing is parsed
        # but then probably nothing moved.
        try:
            del out_dict['atoms']
        except KeyError:
            pass
        try:
            del out_dict['cell']
        except KeyError:
            pass
        try:
            del out_dict['ions_positions_stau']
        except KeyError:
            pass
        try:
            del out_dict['ions_positions_svel']
        except KeyError:
            pass
        try:
            del out_dict['ions_positions_taui']
        except KeyError:
            pass
        # This should not be needed
        try:
            del out_dict['atoms_index_list']
        except KeyError:
            pass
        # This should be already in the input
        try:
            del out_dict['atoms_if_pos_list']
        except KeyError:
            pass
        #
        try:
            del out_dict['ions_positions_force']
        except KeyError:
            pass

        # convert the dictionary into an AiiDA object
        output_params = Dict(dict=out_dict)
        self.out('output_parameters', output_params)

    def get_linkname_trajectory(self):
        """
        Returns the name of the link to the output_structure (None if not present)
        """
        return 'output_trajectory'

    def _generate_sites_ordering(self, raw_species, raw_atoms):
        """
        take the positions of xml and from file.pos of the LAST step and compare them
        """
        # Examples in the comments are for species [Ba, O, Ti]
        # and atoms [Ba, Ti, O, O, O]

        # Dictionary to associate the species name to the idx
        # Example: {'Ba': 1, 'O': 2, 'Ti': 3}
        species_dict = {name: idx for idx, name in zip(raw_species['index'],
                                                       raw_species['type'])}
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
        sorted_indexed_reordering = sorted([(_[1], _[0]) for _ in
                                            enumerate(reordering)])
        reordering_inverse = [_[1] for _ in sorted_indexed_reordering]
        return reordering_inverse

    def _get_reordered_list(self, origlist, reordering):
        """
        Given a list to reorder, a list of integer positions with the new
        order, return the reordered list.
        """
        return [origlist[e] for e in reordering]

    def _get_reordered_array(self, _input, reordering):
        return numpy.array([self._get_reordered_list(i, reordering) for i in _input])
