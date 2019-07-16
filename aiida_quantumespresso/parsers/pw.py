# -*- coding: utf-8 -*-
from __future__ import absolute_import

import numpy
from six.moves import zip

from aiida import orm
from aiida.common import exceptions
from aiida.parsers import Parser

from aiida_quantumespresso.utils.mapping import get_logging_container


class PwParser(Parser):
    """`Parser` implementation for the `PwCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed `PwCalculation` into output nodes.

        Two nodes that are expected are the default 'retrieved' `FolderData` node which will store the retrieved files
        permanently in the repository. The second required node is a filepath under the key `retrieved_temporary_files`
        which should contain the temporary retrieved files.
        """
        dir_with_bands = None
        self.exit_code_xml = None
        self.exit_code_stdout = None
        self.exit_code_parser = None

        try:
            self.retrieved
        except exceptions.NotExistent:
            return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_FOLDER)

        # Look for optional settings input node and potential 'parser_options' dictionary within it
        try:
            parser_options = self.node.inputs.settings.get_dict()[self.get_parser_settings_key()]
        except (KeyError, exceptions.NotExistent):
            parser_options = None

        # Verify that the retrieved_temporary_folder is within the arguments if temporary files were specified
        if self.node.get_attribute('retrieve_temporary_list', None):
            try:
                dir_with_bands = kwargs['retrieved_temporary_folder']
            except KeyError:
                return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_TEMPORARY_FOLDER)

        parameters = self.node.inputs.parameters.get_dict()
        parsed_xml, logs_xml = self.parse_xml(dir_with_bands)
        parsed_stdout, logs_stdout = self.parse_stdout(parameters, parser_options, parsed_xml)

        parsed_bands = parsed_stdout.pop('bands', {})
        parsed_structure = parsed_stdout.pop('structure', {})
        parsed_trajectory = parsed_stdout.pop('trajectory', {})
        parsed_parameters = self.build_output_parameters(parsed_xml, parsed_stdout)

        # Append the last frame of some of the smaller trajectory arrays to the parameters for easy querying
        self.final_trajectory_frame_to_parameters(parsed_parameters, parsed_trajectory)

        # If the parser option 'all_symmetries' is not set to True, we reduce the raw parsed symmetries to safe space
        self.reduce_symmetries(parsed_parameters, parsed_structure, parser_options)

        structure = self.build_output_structure(parsed_structure)
        kpoints = self.build_output_kpoints(parsed_parameters, structure)
        bands = self.build_output_bands(parsed_bands, kpoints)
        array = self.build_output_array(parsed_trajectory)
        trajectory = self.build_output_trajectory(parsed_trajectory, structure)

        # Determine whether the input kpoints were defined as a mesh or as an explicit list
        try:
            self.node.inputs.kpoints.get_kpoints()
        except AttributeError:
            input_kpoints_explicit = False
        else:
            input_kpoints_explicit = True

        # Only attach the `KpointsData` as output if there will be no `BandsData` output and inputs were defined as mesh
        if kpoints and not bands and not input_kpoints_explicit:
            self.out('output_kpoints', kpoints)

        if bands:
            self.out('output_band', bands)

        if trajectory:
            self.out('output_trajectory', trajectory)

        if array:
            self.out('output_array', array)

        if not structure.is_stored:
            self.out('output_structure', structure)

        # Separate the atomic_occupations dictionary in its own node if it is present
        atomic_occupations = parsed_parameters.pop('atomic_occupations', None)
        if atomic_occupations:
            self.out('output_atomic_occupations', orm.Dict(dict=atomic_occupations))

        self.out('output_parameters', orm.Dict(dict=parsed_parameters))

        # Emit the logs returned by the XML and stdout parsing through the logger
        self.emit_logs(logs_stdout, logs_xml)

        # If the both stdout and xml exit codes are set, there was a basic problem with both output files and there
        # is no need to investigate any further.
        if self.exit_code_stdout and self.exit_code_xml:
            return self.exit(self.exit_codes.ERROR_OUTPUT_FILES)

        # First check for specific known problems that can cause a pre-mature termination of the calculation
        exit_code = self.validate_premature_exit(logs_stdout)
        if exit_code:
            return self.exit(exit_code)

        if self.exit_code_stdout:
            return self.exit(self.exit_code_stdout)

        if self.exit_code_xml:
            return self.exit(self.exit_code_xml)

        # First determine issues that can occurr for all calculation types. Note that the generic errors, that are
        # common to all types are done first. If a problem is found there, we return the exit code and don't continue
        for validator in [self.validate_electronic, self.validate_dynamics, self.validate_ionic]:
            exit_code = validator(array, trajectory, parsed_parameters, logs_stdout)
            if exit_code:
                return self.exit(exit_code)

    def validate_premature_exit(self, logs):
        """Analyze problems that will cause a pre-mature termination of the calculation, controlled or not."""
        if 'ERROR_OUT_OF_WALLTIME' in logs['error']:
            return self.exit_codes.ERROR_OUT_OF_WALLTIME

        if 'ERROR_CHARGE_IS_WRONG' in logs['error']:
            return self.exit_codes.ERROR_CHARGE_IS_WRONG

        if 'ERROR_SYMMETRY_NON_ORTHOGONAL_OPERATION' in logs['error']:
            return self.exit_codes.ERROR_SYMMETRY_NON_ORTHOGONAL_OPERATION

    def validate_electronic(self, array, trajectory, parameters, logs):
        """Analyze problems that are specific to `electronic` type calculations: i.e. `scf`, `nscf` and `bands`."""
        if array is None:
            return

        if 'ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED' in logs['error']:
            return self.exit_codes.ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED

    def validate_dynamics(self, array, trajectory, parameters, logs):
        """Analyze problems that are specific to `dynamics` type calculations: i.e. `md` and `vc-md`."""

    def validate_ionic(self, array, trajectory, parameters, logs):
        """Analyze problems that are specific to `ionic` type calculations: i.e. `relax` and `vc-relax`."""
        from aiida_quantumespresso.utils.defaults.calculation import pw
        from aiida_quantumespresso.utils.validation.trajectory import verify_convergence_trajectory

        if trajectory is None:
            return

        electronic_convergence_reached = 'ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED' not in logs.error
        ionic_convergence_reached = 'ERROR_IONIC_CONVERGENCE_NOT_REACHED' not in logs.error
        bfgs_history_failure = 'ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE' in logs.error
        maximum_ionic_steps_reached = 'ERROR_MAXIMUM_IONIC_STEPS_REACHED' in logs.warning
        final_scf = parameters.get('final_scf', False)

        # The electronic self-consistency cycle failed before reaching ionic convergence
        if not ionic_convergence_reached and not electronic_convergence_reached:
            return self.exit_codes.ERROR_IONIC_CYCLE_ELECTRONIC_CONVERGENCE_NOT_REACHED

        # Ionic convergence was not reached because maximum number of steps was exceeded
        if not ionic_convergence_reached and maximum_ionic_steps_reached:
            return self.exit_codes.ERROR_IONIC_CYCLE_EXCEEDED_NSTEP

        # BFGS fails twice in a row in which case QE will print that convergence is reached while it is not
        if ionic_convergence_reached and bfgs_history_failure:

            # If electronic convergence was not reached, this had to have been a `vc-relax` where final SCF failed
            if not electronic_convergence_reached:
                return self.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE

            return self.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE

        threshold_forces = parameters.get('CONTROL', {}).get('forc_conv_thr', pw.forc_conv_thr)
        threshold_stress = parameters.get('CELL', {}).get('press_conv_thr', pw.press_conv_thr)
        thresholds = [threshold_forces, threshold_stress]

        # Here we are converged ionically: if there is no final scf, this was a `relax` run
        if ionic_convergence_reached and not final_scf:

            # Manually verify convergence of total energy and forces
            converged = verify_convergence_trajectory(trajectory, -1, threshold_forces)

            if converged is None or not converged:
                # This should never happen: apparently there was data missing or not converged despite QE saying so
                return self.exit_codes.ERROR_IONIC_CONVERGENCE_NOT_REACHED

        # Ionic convergence was reached in `vc-relax` run, but electronic convergence in final SCF failed
        elif ionic_convergence_reached and not electronic_convergence_reached:

            # Manually verify convergence of total energy and forces
            converged = verify_convergence_trajectory(trajectory, -1, *thresholds)

            if converged is None or not converged:
                # This should never happen: apparently there was data missing or not converged despite QE saying so
                return self.exit_codes.ERROR_IONIC_CONVERGENCE_NOT_REACHED
            else:
                return self.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_FINAL_SCF_FAILED

        # Ionic convergence was reached in `vc-relax` run, and electronic convergence was reached in final SCF
        elif ionic_convergence_reached and electronic_convergence_reached:

            converged_relax = verify_convergence_trajectory(trajectory, -2, *thresholds)
            converged_final = verify_convergence_trajectory(trajectory, -1, *thresholds)

            if converged_final is None or converged_relax is None or not converged_relax:
                # This should never happen: apparently there was data missing or not converged despite QE saying so
                return self.exit_codes.ERROR_IONIC_CONVERGENCE_NOT_REACHED

            if not converged_final:
                # The forces and stresses of ionic cycle are below threshold, but those of the final SCF exceed them.
                # This is not necessarily a problem since the calculation starts from scratch after the variable cell
                # relaxation and the forces and stresses can be slightly different. Still it is useful to distinguish
                # these calculations so we return a special exit code.
                return self.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF

    def exit(self, exit_code):
        """Log the exit message of the give exit code with level `ERROR` and return the exit code.

        This is a utility function if one wants to return from the parse method and automically add the exit message
        associated to the exit code as a log message to the node: e.g. `return self.exit(self.exit_codes.LABEL))`

        :param exit_code: an `ExitCode`
        :return: the exit code
        """
        self.logger.error(exit_code.message)
        return exit_code

    def parse_xml(self, dir_with_bands=None):
        """Parse the XML output file.

        :param dir_with_bands: absolute path to directory containing individual k-point XML files for old XML format.
        :return: tuple of two dictionaries, first with raw parsed data and second with log messages
        """
        from .parse_xml.pw.exceptions import XMLParseError, XMLUnsupportedFormatError
        from .parse_xml.pw.parse import parse_xml

        logs = get_logging_container()
        parsed_data = {}

        object_names = self.retrieved.list_object_names()
        xml_files = [xml_file for xml_file in self.node.process_class.xml_filenames if xml_file in object_names]

        if not xml_files:
            self.exit_code_xml = self.exit_codes.ERROR_OUTPUT_XML_MISSING
            return parsed_data, logs
        elif len(xml_files) > 1:
            self.exit_code_xml = self.exit_codes.ERROR_OUTPUT_XML_MULTIPLE
            return parsed_data, logs

        try:
            with self.retrieved.open(xml_files[0]) as xml_file:
                parsed_data, logs = parse_xml(xml_file, dir_with_bands)
        except IOError:
            self.exit_code_xml = self.exit_codes.ERROR_OUTPUT_XML_READ
        except XMLParseError:
            self.exit_code_xml = self.exit_codes.ERROR_OUTPUT_XML_PARSE
        except XMLUnsupportedFormatError:
            self.exit_code_xml = self.exit_codes.ERROR_OUTPUT_XML_FORMAT
        except Exception:
            import traceback
            traceback.print_exc()
            self.exit_code_xml = self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION

        return parsed_data, logs

    def parse_stdout(self, parameters, parser_options=None, parsed_xml=None):
        """Parse the stdout output file.

        :param parameters: the input parameters dictionary
        :param parser_options: optional dictionary with parser options
        :param parsed_xml: the raw parsed data from the XML output
        :return: tuple of two dictionaries, first with raw parsed data and second with log messages
        """
        from aiida_quantumespresso.parsers.parse_raw.pw import parse_stdout

        logs = get_logging_container()
        parsed_data = {}

        filename_stdout = self.node.get_attribute('output_filename')

        if filename_stdout not in self.retrieved.list_object_names():
            self.exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_MISSING
            return parsed_data, logs

        try:
            stdout = self.retrieved.get_object_content(filename_stdout)
        except IOError:
            self.exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_READ
            return parsed_data, logs

        try:
            parsed_data, logs = parse_stdout(stdout, parameters, parser_options, parsed_xml)
        except Exception:
            import traceback
            traceback.print_exc()
            self.exit_code_stdout = self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION

        # If the stdout was incomplete, most likely the job was interrupted before it could cleanly finish, so the
        # output files are most likely corrupt and cannot be restarted from
        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE' in logs['error']:
            self.exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE

        return parsed_data, logs

    def emit_logs(self, *args):
        """Emit the messages in one or multiple "log dictionaries" through the logger of the parser.

        A log dictionary is expected to have the following structure: each key must correspond to a log level of the
        python logging module, e.g. `error` or `warning` and its values must be a list of string messages. The method
        will loop over all log dictionaries and emit the messages it contains with the log level indicated by the key.

        Example log dictionary structure::

            logs = {
                'warning': ['Could not parse the `etot_threshold` variable from the stdout.'],
                'error': ['Self-consistency was not achieved']
            }

        :param args: log dictionaries
        """
        ignore = [
            'Error while parsing ethr.',
            'DEPRECATED: symmetry with ibrav=0, use correct ibrav instead'
        ]

        for logs in args:
            for level, messages in logs.items():
                for message in messages:

                    if message is None:
                        continue

                    stripped = message.strip()

                    if not stripped or stripped in ignore:
                        continue

                    try:
                        getattr(self.logger, level)(stripped)
                    except AttributeError:
                        pass

    def build_output_parameters(self, parsed_stdout, parsed_xml):
        """Build the dictionary of output parameters from the raw parsed data.

        The output parameters are based on the union of raw parsed data from the XML and stdout output files.
        Currently, if both raw parsed data dictionaries contain the same key, the stdout version takes precedence, but
        this should not occur as the `parse_stdout` method should already have solved these conflicts.

        :param parsed_stdout: the raw parsed data dictionary from the stdout output file
        :param parsed_xml: the raw parsed data dictionary from the XML output file
        :return: the union of the two parsed raw and information about the parser
        """
        from aiida_quantumespresso.parsers import get_parser_info

        parsed_info = get_parser_info(parser_info_template='aiida-quantumespresso parser pw.x v{}')

        for key in list(parsed_stdout.keys()):
            if key in list(parsed_xml.keys()):
                if key == 'fermi_energy' or key == 'fermi_energy_units':  # an exception for the (only?) key that may be found on both
                    del parsed_stdout[key]
                else:
                    raise AssertionError('{} found in both dictionaries, values: {} vs. {}'.format(
                        key, parsed_stdout[key], parsed_xml[key]))  # this shouldn't happen!

        parameters = dict(list(parsed_xml.items()) + list(parsed_stdout.items()) + list(parsed_info.items()))

        return parameters

    def build_output_structure(self, parsed_structure):
        """Build the output structure from the raw parsed data.

        :param parsed_structure: the dictionary with raw parsed structure data
        :return: a new `StructureData` created from the parsed data iff the calculation type produces a new structure
            and the parsed data contained a cell definition. In all other cases, the input structure will be returned.
        """
        from aiida_quantumespresso.parsers import convert_qe2aiida_structure

        type_calc = self.node.inputs.parameters.get_dict()['CONTROL']['calculation']

        if type_calc not in ['relax', 'vc-relax', 'md', 'vc-md'] or 'cell' not in list(parsed_structure.keys()):
            return self.node.inputs.structure

        return convert_qe2aiida_structure(parsed_structure, self.node.inputs.structure)

    def build_output_array(self, parsed_trajectory):
        """Build the output array from the raw parsed trajectory data.

        :param parsed_trajectory: the raw parsed trajectory data
        :return: an `ArrayData` node
        """
        if not parsed_trajectory or 'atomic_positions_relax' in parsed_trajectory:
            return None

        array = orm.ArrayData()

        for frame in parsed_trajectory.items():
            array.set_array(frame[0], numpy.array(frame[1]))

        return array

    def build_output_trajectory(self, parsed_trajectory, structure):
        """Build the output trajectory from the raw parsed trajectory data.

        :param parsed_trajectory: the raw parsed trajectory data
        :return: a `TrajectoryData` or None
        """
        if not parsed_trajectory or 'atomic_positions_relax' not in parsed_trajectory:
            return None

        positions = numpy.array(parsed_trajectory.pop('atomic_positions_relax'))
        try:
            cells = numpy.array(parsed_trajectory.pop('lattice_vectors_relax'))
            # if the cell is only printed once, the MD/relax was at fixed cell
            if len(cells) == 1 and len(positions) > 1:
                cells = numpy.array([cells[0]] * len(positions))
        except KeyError:
            # The cell is never printed, the MD/relax was at fixed cell
            cells = numpy.array([structure.cell] * len(positions))

        symbols = [str(site.kind_name) for site in structure.sites]
        stepids = numpy.arange(len(positions))

        trajectory = orm.TrajectoryData()
        trajectory.set_trajectory(
            stepids=stepids,
            cells=cells,
            symbols=symbols,
            positions=positions,
        )

        for frame in parsed_trajectory.items():
            trajectory.set_array(frame[0], numpy.array(frame[1]))

        return trajectory

    def build_output_kpoints(self, parsed_parameters, structure):
        """Build the output kpoints from the raw parsed data.

        :param parsed_parameters: the raw parsed data
        :return: a `KpointsData` or None
        """
        k_points_list = parsed_parameters.pop('k_points', None)
        k_points_units = parsed_parameters.pop('k_points_units', None)
        k_points_weights_list = parsed_parameters.pop('k_points_weights', None)

        if k_points_list is None or k_points_weights_list is None:
            return None

        if k_points_units != '1 / angstrom':
            self.logger.error('Error in kpoints units (should be cartesian)')
            self.exit_code_parser.ERROR_INVALID_KPOINT_UNITS

        kpoints = orm.KpointsData()
        kpoints.set_cell_from_structure(structure)
        kpoints.set_kpoints(k_points_list, cartesian=True, weights=k_points_weights_list)

        return kpoints

    def build_output_bands(self, parsed_bands, parsed_kpoints=None):
        """Build the output bands from the raw parsed bands data.

        :param parsed_bands: the raw parsed bands data
        :param parsed_kpoints: the `KpointsData` to use for the bands
        :return: a `BandsData` or None
        """
        if not parsed_bands or not parsed_kpoints:
            return

        # Correct the occupation for nspin=1 calculations where Quantum ESPRESSO populates each band only halfway
        if len(parsed_bands['occupations']) > 1:
            occupations = parsed_bands['occupations']
        else:
            occupations = 2. * numpy.array(parsed_bands['occupations'][0])

        if len(parsed_bands['bands']) > 1:
            bands_energies = parsed_bands['bands']
        else:
            bands_energies = parsed_bands['bands'][0]

        bands = orm.BandsData()
        bands.set_kpointsdata(parsed_kpoints)
        bands.set_bands(bands_energies, units=parsed_bands['bands_units'], occupations=occupations)

        return bands

    def get_parser_settings_key(self):
        """
        Return the name of the key to be used in the calculation settings, that
        contains the dictionary with the parser_options
        """
        return 'parser_options'

    def final_trajectory_frame_to_parameters(self, parameters, parsed_trajectory):
        """Move the last frame of certain properties from the `TrajectoryData` to the outputs parameters.

        This makes these properties queryable.
        """
        ignore_keys = [
            'atomic_charges',
            'atomic_magnetic_moments',
            'atomic_positions_relax',
            'atomic_species_name',
            'forces',
            'stress',
            'lattice_vectors_relax',
        ]

        # Have to iterate over a static list of keys since we are mutating `parsed_trajectory` within the loop
        for property_key in list(parsed_trajectory.keys()):

            property_values = parsed_trajectory[property_key]

            if property_key in ignore_keys:
                continue

            parameters[property_key] = property_values[-1]

            # Delete properties whose values list has but a single value
            if len(property_values) == 1:
                parsed_trajectory.pop(property_key)

    def reduce_symmetries(self, parsed_parameters, parsed_structure, parser_options=None):
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
        from aiida_quantumespresso.utils.linalg import are_matrices_equal

        all_symmetries = False if parser_options is None else parser_options.get('all_symmetries', False)

        if all_symmetries or 'cell' not in parsed_structure:
            return

        cell = parsed_structure['cell']['lattice_vectors']
        cell_T = numpy.transpose(cell)
        cell_Tinv = numpy.linalg.inv(cell_T)
        possible_symmetries = self._get_qe_symmetry_list()

        # for symmetry_type in ['symmetries', 'lattice_symmetries']:  # crystal vs. lattice symmetries
        for symmetry_type in ['symmetries']:  # crystal vs. lattice symmetries
            if symmetry_type in list(parsed_parameters.keys()):
                try:
                    old_symmetries = parsed_parameters[symmetry_type]
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

                    parsed_parameters[symmetry_type] = new_symmetries  # and overwrite the old one
                except KeyError:
                    self.logger.warning("key '{}' is not present in raw output dictionary".format(symmetry_type))
            else:
                # backwards-compatiblity: 'lattice_symmetries' is not created in older versions of the parser
                if symmetry_type != 'lattice_symmetries':
                    self.logger.warning("key '{}' is not present in raw output dictionary".format(symmetry_type))

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
