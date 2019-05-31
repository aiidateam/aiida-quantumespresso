# -*- coding: utf-8 -*-
from __future__ import absolute_import

import copy
import numpy
import six
from six.moves import zip

from aiida import orm
from aiida.common import exceptions
from aiida.engine import ExitCode
from aiida.parsers import Parser
from aiida_quantumespresso.parsers import convert_qe2aiida_structure
from aiida_quantumespresso.parsers.parse_raw.pw import parse_raw_output
from aiida_quantumespresso.utils.linalg import are_matrices_equal


class PwParser(Parser):
    """`Parser` implementation for the `PwCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the output nodes for a PwCalculations from a dictionary of retrieved nodes.

        Two nodes that are expected are the default 'retrieved' `FolderData` node which will store the retrieved files
        permanently in the repository. The second required node is a filepath under the key `retrieved_temporary_files`
        which should contain the temporary retrieved files.
        """
        try:
            retrieved = self.retrieved
        except exceptions.NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        input_parameters = self.node.inputs.parameters.get_dict()

        # Look for optional settings input node and potential 'parser_options' dictionary within it
        try:
            settings = self.node.inputs.settings.get_dict()
            parser_options = settings[self.get_parser_settings_key()]
        except (AttributeError, KeyError, exceptions.NotExistent):
            settings = {}
            parser_options = {}

        # Verify that the retrieved_temporary_folder is within the arguments if temporary files were specified
        if self.node.get_attribute('retrieve_temporary_list', None):
            try:
                temporary_folder = kwargs['retrieved_temporary_folder']
                dir_with_bands = temporary_folder
            except KeyError:
                self.logger.error('the `retrieved_temporary_folder` was not passed as an argument')
                return self.exit_codes.ERROR_NO_RETRIEVED_TEMPORARY_FOLDER
        else:
            dir_with_bands = None

        list_of_files = retrieved.list_object_names()

        try:
            filename_stdout = self.node.get_attribute('output_filename')
            stdout = retrieved.get_object_content(filename_stdout)
        except IOError:
            self.logger.error("The required standard output file '{}' could not be read".format(filename_stdout))
            return self.exit_codes.ERROR_READING_OUTPUT_FILE

        # The xml file is required for successful parsing, but if it does not exist, we still try to parse the out file
        xml_files = [xml_file for xml_file in self.node.process_class.xml_filenames if xml_file in list_of_files]

        if not xml_files:
            self.logger.error('no XML output files found, which is required for parsing')
            return self.exit_codes.ERROR_MISSING_XML_FILE
        elif len(xml_files) > 1:
            self.logger.error('more than one XML file retrieved, which should never happen')
            return self.exit_codes.ERROR_MULTIPLE_XML_FILES

        with retrieved.open(xml_files[0]) as xml_file:
            parsing_args = [stdout, xml_file, input_parameters, parser_options, dir_with_bands]
            try:
                parsed_data, logs = parse_raw_output(*parsing_args)
            except Exception:
                import traceback
                traceback.print_exc()
                self.logger.error('parsing of output files raised an exception')
                return self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION

        # Log the messages returned by the parser through the logger of the parser. This will ensure that they get
        # attached to the report of the corresponding node.
        for level, messages in logs.items():
            for message in messages:
                getattr(self.logger, level)(message)

        bands_data = parsed_data['bands']
        parameters = parsed_data['parameters']
        structure_data = parsed_data['structure']
        trajectory_data = parsed_data['trajectory']

        # I add in the parsed_stdout all the last elements of trajectory_data values,
        # except for some large arrays, that I will likely never query.
        skip_keys = ['forces', 'atomic_magnetic_moments', 'atomic_charges', 'lattice_vectors_relax',
                     'atomic_positions_relax', 'atomic_species_name']
        tmp_trajectory_data = copy.copy(trajectory_data)
        for x in six.iteritems(tmp_trajectory_data):
            if x[0] in skip_keys:
                continue
            parameters[x[0]] = x[1][-1]
            if len(x[1]) == 1:  # delete eventual keys that are not arrays (scf cycles)
                trajectory_data.pop(x[0])
            # note: if an array is empty, there will be KeyError

        # As the k points are an array that is rather large, and again it's not something I'm going to parse likely
        # since it's an info mainly contained in the input file, I move it to the trajectory data
        for key in ['k_points', 'k_points_weights']:
            try:
                trajectory_data[key] = parameters.pop(key)
            except KeyError:
                pass

        # If the parser option 'all_symmetries' is not set to True, we reduce the raw parsed symmetries to safe space
        all_symmetries = parser_options.get('all_symmetries', False)

        if not all_symmetries and 'cell' in structure_data:
            self.reduce_symmetries(parameters, structure_data)

        # I eventually save the new structure. structure_data is unnecessary after this
        input_structure = self.node.get_incoming(link_label_filter='structure').one().node
        type_calc = input_parameters['CONTROL']['calculation']
        output_structure = input_structure
        if type_calc in ['relax', 'vc-relax', 'md', 'vc-md']:
            if 'cell' in list(structure_data.keys()):
                output_structure = convert_qe2aiida_structure(structure_data, input_structure=input_structure)
                self.out('output_structure', output_structure)

        k_points_list = trajectory_data.pop('k_points', None)
        k_points_weights_list = trajectory_data.pop('k_points_weights', None)

        if k_points_list is not None:

            if parameters.pop('k_points_units') != '1 / angstrom':
                self.logger.error('Error in kpoints units (should be cartesian)')
                return self.exit_codes.ERROR_INVALID_KPOINT_UNITS

            kpoints_from_output = orm.KpointsData()
            kpoints_from_output.set_cell_from_structure(output_structure)
            kpoints_from_output.set_kpoints(k_points_list, cartesian=True, weights=k_points_weights_list)
            kpoints_from_input = self.node.inputs.kpoints

            if not bands_data:
                try:
                    kpoints_from_input.get_kpoints()
                except AttributeError:
                    self.out('output_kpoints', kpoints_from_output)

            # Converting bands into a BandsData object (including the kpoints)
            if bands_data:
                kpoints_for_bands = kpoints_from_output

                try:
                    kpoints_from_input.get_kpoints()
                    kpoints_for_bands.labels = kpoints_from_input.labels
                except (AttributeError, ValueError, TypeError):
                    # AttributeError: no list of kpoints in input
                    # ValueError: labels from input do not match the output
                    #      list of kpoints (some kpoints are missing)
                    # TypeError: labels are not set, so kpoints_from_input.labels=None
                    pass

                # Get the bands occupations and correct the occupations of QE:
                # If it computes only one component, it occupies it with half number of electrons
                # TODO: is this something we really want to correct?
                # Probably yes for backwards compatiblity, but it might not be intuitive.
                if len(bands_data['occupations']) > 1:
                    the_occupations = bands_data['occupations']
                else:
                    the_occupations = 2. * numpy.array(bands_data['occupations'][0])

                if len(bands_data['bands']) > 1:
                    bands_energies = bands_data['bands']
                else:
                    bands_energies = bands_data['bands'][0]

                # TODO: replicate this in the new parser!
                the_bands_data = orm.BandsData()
                the_bands_data.set_kpointsdata(kpoints_for_bands)
                the_bands_data.set_bands(bands_energies,
                                         units=bands_data['bands_units'],
                                         occupations=the_occupations)

                self.out('output_band', the_bands_data)
                parameters['linknames_band'] = ['output_band']
                # TODO: replicate this in the new parser!

        # Separate the atomic_occupations dictionary in its own node if it is present
        atomic_occupations = parameters.pop('atomic_occupations', None)
        if atomic_occupations:
            self.out('output_atomic_occupations', orm.Dict(dict=atomic_occupations))

        self.out('output_parameters', orm.Dict(dict=parameters))

        if trajectory_data:
            try:
                positions = numpy.array(trajectory_data.pop('atomic_positions_relax'))
                try:
                    cells = numpy.array(trajectory_data.pop('lattice_vectors_relax'))
                    # if the cell is only printed once, the MD/relax was at fixed cell
                    if len(cells) == 1 and len(positions) > 1:
                        cells = numpy.array([cells[0]] * len(positions))
                # if KeyError (the cell is never printed), the MD/relax was at fixed cell
                except KeyError:
                    cells = numpy.array([input_structure.cell] * len(positions))

                symbols = [str(i.kind_name) for i in input_structure.sites]
                stepids = numpy.arange(len(positions))  # a growing integer per step
                # I will insert time parsing when they fix their issues about time
                # printing (logic is broken if restart is on)

                traj = orm.TrajectoryData()
                traj.set_trajectory(
                    stepids=stepids,
                    cells=cells,
                    symbols=symbols,
                    positions=positions,
                )

                for x in six.iteritems(trajectory_data):
                    traj.set_array(x[0], numpy.array(x[1]))
                self.out('output_trajectory', traj)

            except KeyError:
                # forces, atomic charges and atomic mag. moments, in scf calculation (when outputed)
                arraydata = orm.ArrayData()
                for x in six.iteritems(trajectory_data):
                    arraydata.set_array(x[0], numpy.array(x[1]))
                self.out('output_array', arraydata)

        return ExitCode(0)

    def get_parser_settings_key(self):
        """
        Return the name of the key to be used in the calculation settings, that
        contains the dictionary with the parser_options
        """
        return 'parser_options'

    def get_extended_symmetries(self):
        """Return the extended dictionary of symmetries based on reduced symmetries stored in output parameters."""
        possible_symmetries = self._get_qe_symmetry_list()
        parameters = self.node.get_outgoing(node_class=orm.Dict).get_node_by_label('output_parameters')

        symmetries_extended = []
        symmetries_reduced = parameters.get_dict()['symmetries']  # rimetti lo zero

        for element in symmetries_reduced:

            symmetry = {}

            for keys in ['t_rev', 'equivalent_ions', 'fractional_translation']:
                try:
                    symmetry[keys] = element[keys]
                except KeyError:
                    pass

            # expand the rest
            symmetry['name'] = possible_symmetries[element['symmetry_number']]['name']
            symmetry['rotation'] = possible_symmetries[element['symmetry_number']]['matrix']
            symmetry['inversion'] = possible_symmetries[element['symmetry_number']]['inversion']

            symmetries_extended.append(symmetry)

        return symmetries_extended

    def reduce_symmetries(self, out_dict, structure_data):
        """Reduce the symmetry information parsed from the output to save space.

        In the standard output, each symmetry operation print two rotation matrices:

            * S_cryst^T: matrix in crystal coordinates, transposed
            * S_cart: matrix in cartesian coordinates,

        The XML files only print one matrix:

            * S_cryst: matrix in crystal coordinates

        The raw parsed symmetry information from the XML is large and will load the database heavily if stored as
        is for each calculation. Instead, we will map these dictionaries onto a static dictionary of rotation
        matrices generated by the _get_qe_symmetry_list static method. This dictionary will return the rotation
        matrices in cartesian coordinates, i.e. S_cart. In order to compare the raw matrices from the XML to these
        static matrices we have to convert S_cryst into S_cart. We derive here how that is done:

            S_cryst * v_cryst = v_cryst'

        where v_cryst' is the rotated vector v_cryst under S_cryst
        We define `cell` where cell vectors are rows. Converting a vector from crystal to cartesian
        coordinates is defined as:

            cell^T * v_cryst = v_cart

        The inverse of this operation is defined as

            v_cryst = cell^Tinv * v_cart

        Replacing the last equation into the first we find:

            S_cryst * cell^Tinv * v_cart = cell^Tinv * v_cart'

        Multiply on the left with cell^T gives:

            cell^T * S_cryst * cell^Tinv * v_cart = v_cart'

        which can be rewritten as:

            S_cart * v_cart = v_cart'

        where:

            S_cart = cell^T * S_cryst * cell^Tinv

        We compute here the transpose and its inverse of the structure cell basis, which is needed to transform
        the parsed rotation matrices, which are in crystal coordinates, to cartesian coordinates, which are the
        matrices that are returned by the _get_qe_symmetry_list staticmethod
        """
        cell = structure_data['cell']['lattice_vectors']
        cell_T = numpy.transpose(cell)
        cell_Tinv = numpy.linalg.inv(cell_T)
        possible_symmetries = self._get_qe_symmetry_list()

        for symmetry_type in ['symmetries', 'lattice_symmetries']:  # crystal vs. lattice symmetries
            if symmetry_type in list(out_dict.keys()):
                try:
                    old_symmetries = out_dict[symmetry_type]
                    new_symmetries = []
                    for this_sym in old_symmetries:
                        name = this_sym['name'].strip()
                        for i, this in enumerate(possible_symmetries):
                            # Since we do an exact comparison we strip the string name from whitespace
                            # and as soon as it is matched, we break to prevent it from matching another
                            if name == this['name'].strip():
                                index = i
                                break
                        else:
                            index = None
                            self.logger.error('Symmetry {} not found'.format(name))

                        new_dict = {}
                        if index is not None:
                            # The raw parsed rotation matrix is in crystal coordinates, whereas the mapped rotation
                            # in possible_symmetries is in cartesian coordinates. To allow them to be compared
                            # to make sure we matched the correct rotation symmetry, we first convert the parsed matrix
                            # to cartesian coordinates. For explanation of the method, see comment above.
                            rotation_cryst = this_sym['rotation']
                            rotation_cart_new = possible_symmetries[index]['matrix']
                            rotation_cart_old = numpy.dot(cell_T, numpy.dot(rotation_cryst, cell_Tinv))

                            inversion = possible_symmetries[index]['inversion']
                            if not are_matrices_equal(rotation_cart_old, rotation_cart_new, swap_sign_matrix_b=inversion):
                                self.logger.error('Mapped rotation matrix {} does not match the original rotation {}'
                                    .format(rotation_cart_new, rotation_cart_old))
                                new_dict['all_symmetries'] = this_sym
                            else:
                                # Note: here I lose the information about equivalent ions and fractional_translation
                                # since I don't copy them to new_dict (but they can be reconstructed).
                                new_dict['t_rev'] = this_sym['t_rev']
                                new_dict['symmetry_number'] = index
                        else:
                            new_dict['all_symmetries'] = this_sym

                        new_symmetries.append(new_dict)

                    out_dict[symmetry_type] = new_symmetries  # and overwrite the old one
                except KeyError:
                    self.logger.warning("key '{}' is not present in raw output dictionary".format(symmetry_type))
            else:
                # backwards-compatiblity: 'lattice_symmetries' is not created in older versions of the parser
                if symmetry_type != 'lattice_symmetries':
                    self.logger.warning("key '{}' is not present in raw output dictionary".format(symmetry_type))

    @staticmethod
    def _get_qe_symmetry_list():
        """
        Hard coded names and rotation matrices + inversion from QE v 5.0.2
        Function for Parser class usage only.

        :return: a list of dictionaries, each containing name (string),
            inversion (boolean) and matrix (list of lists)
        """
        sin3 = 0.866025403784438597
        cos3 = 0.5
        msin3 = -0.866025403784438597
        mcos3 = -0.5

        # 32 rotations that are checked + inversion taken from symm_base.f90 from the QE source code
        # They are in Fortran format and therefore transposed with respect to the default python format
        transposed_matrices_cartesian = [
            [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            [[-1., 0., 0.], [0., -1., 0.], [0., 0., 1.]],
            [[-1., 0., 0.], [0., 1., 0.], [0., 0., -1.]],
            [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
            [[0., 1., 0.], [1., 0., 0.], [0., 0., -1.]],
            [[0., -1., 0.], [-1., 0., 0.], [0., 0., -1.]],
            [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
            [[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]],
            [[0., 0., 1.], [0., -1., 0.], [1., 0., 0.]],
            [[0., 0., -1.], [0., -1., 0.], [-1., 0., 0.]],
            [[0., 0., -1.], [0., 1., 0.], [1., 0., 0.]],
            [[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]],
            [[-1., 0., 0.], [0., 0., 1.], [0., 1., 0.]],
            [[-1., 0., 0.], [0., 0., -1.], [0., -1., 0.]],
            [[1., 0., 0.], [0., 0., -1.], [0., 1., 0.]],
            [[1., 0., 0.], [0., 0., 1.], [0., -1., 0.]],
            [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
            [[0., 0., -1.], [-1., 0., 0.], [0., 1., 0.]],
            [[0., 0., -1.], [1., 0., 0.], [0., -1., 0.]],
            [[0., 0., 1.], [-1., 0., 0.], [0., -1., 0.]],
            [[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]],
            [[0., -1., 0.], [0., 0., -1.], [1., 0., 0.]],
            [[0., -1., 0.], [0., 0., 1.], [-1., 0., 0.]],
            [[0., 1., 0.], [0., 0., -1.], [-1., 0., 0.]],
            [[cos3, sin3, 0.], [msin3, cos3, 0.], [0., 0., 1.]],
            [[cos3, msin3, 0.], [sin3, cos3, 0.], [0., 0., 1.]],
            [[mcos3, sin3, 0.], [msin3, mcos3, 0.], [0., 0., 1.]],
            [[mcos3, msin3, 0.], [sin3, mcos3, 0.], [0., 0., 1.]],
            [[cos3, msin3, 0.], [msin3, mcos3, 0.], [0., 0., -1.]],
            [[cos3, sin3, 0.], [sin3, mcos3, 0.], [0., 0., -1.]],
            [[mcos3, msin3, 0.], [msin3, cos3, 0.], [0., 0., -1.]],
            [[mcos3, sin3, 0.], [sin3, cos3, 0.], [0., 0., -1.]],
        ]

        # Names for the 32 matrices, with and without inversion
        matrices_name = [
            'identity                                     ',
            '180 deg rotation - cart. axis [0,0,1]        ',
            '180 deg rotation - cart. axis [0,1,0]        ',
            '180 deg rotation - cart. axis [1,0,0]        ',
            '180 deg rotation - cart. axis [1,1,0]        ',
            '180 deg rotation - cart. axis [1,-1,0]       ',
            ' 90 deg rotation - cart. axis [0,0,-1]       ',
            ' 90 deg rotation - cart. axis [0,0,1]        ',
            '180 deg rotation - cart. axis [1,0,1]        ',
            '180 deg rotation - cart. axis [-1,0,1]       ',
            ' 90 deg rotation - cart. axis [0,1,0]        ',
            ' 90 deg rotation - cart. axis [0,-1,0]       ',
            '180 deg rotation - cart. axis [0,1,1]        ',
            '180 deg rotation - cart. axis [0,1,-1]       ',
            ' 90 deg rotation - cart. axis [-1,0,0]       ',
            ' 90 deg rotation - cart. axis [1,0,0]        ',
            '120 deg rotation - cart. axis [-1,-1,-1]     ',
            '120 deg rotation - cart. axis [-1,1,1]       ',
            '120 deg rotation - cart. axis [1,1,-1]       ',
            '120 deg rotation - cart. axis [1,-1,1]       ',
            '120 deg rotation - cart. axis [1,1,1]        ',
            '120 deg rotation - cart. axis [-1,1,-1]      ',
            '120 deg rotation - cart. axis [1,-1,-1]      ',
            '120 deg rotation - cart. axis [-1,-1,1]      ',
            ' 60 deg rotation - cryst. axis [0,0,1]       ',
            ' 60 deg rotation - cryst. axis [0,0,-1]      ',
            '120 deg rotation - cryst. axis [0,0,1]       ',
            '120 deg rotation - cryst. axis [0,0,-1]      ',
            '180 deg rotation - cryst. axis [1,-1,0]      ',
            '180 deg rotation - cryst. axis [2,1,0]       ',
            '180 deg rotation - cryst. axis [0,1,0]       ',
            '180 deg rotation - cryst. axis [1,1,0]       ',
            'inversion                                    ',
            'inv. 180 deg rotation - cart. axis [0,0,1]   ',
            'inv. 180 deg rotation - cart. axis [0,1,0]   ',
            'inv. 180 deg rotation - cart. axis [1,0,0]   ',
            'inv. 180 deg rotation - cart. axis [1,1,0]   ',
            'inv. 180 deg rotation - cart. axis [1,-1,0]  ',
            'inv.  90 deg rotation - cart. axis [0,0,-1]  ',
            'inv.  90 deg rotation - cart. axis [0,0,1]   ',
            'inv. 180 deg rotation - cart. axis [1,0,1]   ',
            'inv. 180 deg rotation - cart. axis [-1,0,1]  ',
            'inv.  90 deg rotation - cart. axis [0,1,0]   ',
            'inv.  90 deg rotation - cart. axis [0,-1,0]  ',
            'inv. 180 deg rotation - cart. axis [0,1,1]   ',
            'inv. 180 deg rotation - cart. axis [0,1,-1]  ',
            'inv.  90 deg rotation - cart. axis [-1,0,0]  ',
            'inv.  90 deg rotation - cart. axis [1,0,0]   ',
            'inv. 120 deg rotation - cart. axis [-1,-1,-1]',
            'inv. 120 deg rotation - cart. axis [-1,1,1]  ',
            'inv. 120 deg rotation - cart. axis [1,1,-1]  ',
            'inv. 120 deg rotation - cart. axis [1,-1,1]  ',
            'inv. 120 deg rotation - cart. axis [1,1,1]   ',
            'inv. 120 deg rotation - cart. axis [-1,1,-1] ',
            'inv. 120 deg rotation - cart. axis [1,-1,-1] ',
            'inv. 120 deg rotation - cart. axis [-1,-1,1] ',
            'inv.  60 deg rotation - cryst. axis [0,0,1]  ',
            'inv.  60 deg rotation - cryst. axis [0,0,-1] ',
            'inv. 120 deg rotation - cryst. axis [0,0,1]  ',
            'inv. 120 deg rotation - cryst. axis [0,0,-1] ',
            'inv. 180 deg rotation - cryst. axis [1,-1,0] ',
            'inv. 180 deg rotation - cryst. axis [2,1,0]  ',
            'inv. 180 deg rotation - cryst. axis [0,1,0]  ',
            'inv. 180 deg rotation - cryst. axis [1,1,0]  '
        ]

        rotations = []

        for key, value in zip(matrices_name[:len(transposed_matrices_cartesian)], transposed_matrices_cartesian):
            rotations.append({'name': key, 'matrix': numpy.transpose(value), 'inversion': False})

        for key, value in zip(matrices_name[len(transposed_matrices_cartesian):], transposed_matrices_cartesian):
            rotations.append({'name': key, 'matrix': numpy.transpose(value), 'inversion': True})

        return rotations
