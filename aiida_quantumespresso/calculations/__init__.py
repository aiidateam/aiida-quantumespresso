# -*- coding: utf-8 -*-
from __future__ import absolute_import

import abc
import io
import os

from aiida import orm
from aiida.common import datastructures, exceptions
from aiida.common.lang import classproperty
from aiida.engine import CalcJob
from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry
import six
from six.moves import zip


class BasePwCpInputGenerator(CalcJob):

    _PSEUDO_SUBFOLDER = './pseudo/'
    _OUTPUT_SUBFOLDER = './out/'
    _PREFIX = 'aiida'
    _INPUT_FILE_NAME = 'aiida.in'
    _OUTPUT_FILE_NAME = 'aiida.out'
    _DATAFILE_XML_PRE_6_2 = 'data-file.xml'
    _DATAFILE_XML_POST_6_2 = 'data-file-schema.xml'
    _ENVIRON_INPUT_FILE_NAME = 'environ.in'

    # Additional files that should always be retrieved for the specific plugin
    _internal_retrieve_list = []

    # Name lists to print by calculation type
    _automatic_namelists = {}

    # Blocked keywords that are to be specified in the subclass
    _blocked_keywords = {}

    # In restarts, will not copy but use symlinks
    _default_symlink_usage = True

    # In restarts, it will copy from the parent the following
    _restart_copy_from = os.path.join(_OUTPUT_SUBFOLDER, '*')

    # In restarts, it will copy the previous folder in the following one
    _restart_copy_to = _OUTPUT_SUBFOLDER

    # Default verbosity; change in subclasses
    _default_verbosity = 'high'

    _use_kpoints = False

    @classproperty
    def xml_filenames(cls):
        """Return a list of XML output filenames that can be written by a calculation.

        Note that this includes all potential filenames across all known versions of Quantum ESPRESSO
        """
        return [cls._DATAFILE_XML_POST_6_2, cls._DATAFILE_XML_PRE_6_2]

    @abc.abstractmethod
    @classproperty
    def xml_filepaths(cls):
        """Return a list of XML output filepaths relative to the remote working directory that should be retrieved."""
        pass

    @classmethod
    def define(cls, spec):
        super(BasePwCpInputGenerator, cls).define(spec)
        spec.input('metadata.options.withmpi', valid_type=bool, default=True)  # Override default withmpi=False
        spec.input('structure', valid_type=orm.StructureData, help='')
        spec.input('parameters', valid_type=orm.Dict, help='')
        spec.input('settings', valid_type=orm.Dict, required=False, help='')
        spec.input('parent_folder', valid_type=orm.RemoteData, required=False, help='')
        spec.input('vdw_table', valid_type=orm.SinglefileData, required=False, help='')
        spec.input_namespace('pseudos', valid_type=orm.UpfData, dynamic=True, help='')

    def prepare_for_submission(self, folder):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        if 'settings' in self.inputs:
            settings = _uppercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
        else:
            settings = {}

        # Check that a pseudo potential was specified for each kind present in the `StructureData`
        kinds = [kind.name for kind in self.inputs.structure.kinds]
        if set(kinds) != set(self.inputs.pseudos.keys()):
            raise exceptions.InputValidationError(
                'Mismatch between the defined pseudos and the list of kinds of the structure.\n'
                'Pseudos: {};\nKinds: {}'.format(', '.join(list(self.inputs.pseudos.keys())), ', '.join(list(kinds))))

        local_copy_list = []
        remote_copy_list = []
        remote_symlink_list = []

        # Create the subfolder that will contain the pseudopotentials
        folder.get_subfolder(self._PSEUDO_SUBFOLDER, create=True)
        # Create the subfolder for the output data (sometimes Quantum ESPRESSO codes crash if the folder does not exist)
        folder.get_subfolder(self._OUTPUT_SUBFOLDER, create=True)

        # If present, add also the Van der Waals table to the pseudo dir. Note that the name of the table is not checked
        # but should be the one expected by Quantum ESPRESSO.
        if 'vdw_table' in self.inputs:
            src_path = self.inputs.vdw_table.get_file_abs_path()
            dst_path = os.path.join(self._PSEUDO_SUBFOLDER, os.path.split(self.inputs.vdw_table.get_file_abs_path())[1])
            local_copy_list.append((src_path, dst_path))

        if 'hubbard_file' in self.inputs:
            src_path = self.inputs.hubbard_file.get_file_abs_path()
            dst_path = self.input_file_name_hubbard_file
            local_copy_list.append((src_path, dst_path))

        arguments = [
            self.inputs.parameters,
            settings,
            self.inputs.pseudos,
            self.inputs.structure,
        ]
        if self._use_kpoints:
            arguments.append(self.inputs.kpoints)
        input_filecontent, local_copy_pseudo_list = self._generate_PWCPinputdata(*arguments)
        local_copy_list += local_copy_pseudo_list

        input_filename = folder.get_abs_path(self._INPUT_FILE_NAME)
        with io.open(input_filename, 'w') as handle:
            handle.write(input_filecontent)

        # operations for restart
        symlink = settings.pop('PARENT_FOLDER_SYMLINK', self._default_symlink_usage)  # a boolean
        if symlink:
            if 'parent_folder' in self.inputs:
                # I put the symlink to the old parent ./out folder
                remote_symlink_list.append((
                    self.inputs.parent_folder.computer.uuid,
                    os.path.join(self.inputs.parent_folder.get_remote_path(), self._restart_copy_from),
                    self._restart_copy_to
                ))
        else:
            # copy remote output dir, if specified
            if 'parent_folder' in self.inputs:
                remote_copy_list.append((
                    self.inputs.parent_folder.computer.uuid,
                    os.path.join(self.inputs.parent_folder.get_remote_path(), self._restart_copy_from),
                    self._restart_copy_to
                ))

        # Create an `.EXIT` file if `only_initialization` flag in `settings` is set to `True`
        if settings.pop('ONLY_INITIALIZATION', False):
            exit_filename = folder.get_abs_path('{}.EXIT'.format(self._PREFIX))
            with open(exit_filename, 'w') as handle:
                handle.write('\n')

        # Check if specific inputs for the ENVIRON module where specified
        environ_namelist = settings.pop('ENVIRON', None)
        if environ_namelist is not None:
            if not isinstance(environ_namelist, dict):
                raise exceptions.InputValidationError("ENVIRON namelist should be specified as a dictionary")
            # We first add the environ flag to the command-line options (if not already present)
            try:
                if '-environ' not in settings['CMDLINE']:
                    settings['CMDLINE'].append('-environ')
            except KeyError:
                settings['CMDLINE'] = ['-environ']
            # To create a mapping from the species to an incremental fortran 1-based index
            # we use the alphabetical order as in the inputdata generation
            kind_names = sorted([kind.name for kind in self.inputs.structure.kinds])
            mapping_species = {kind_name: (index + 1) for index, kind_name in enumerate(kind_names)}
            environ_input_filename = folder.get_abs_path(self._ENVIRON_INPUT_FILE_NAME)

            with open(environ_input_filename, 'w') as handle:
                handle.write('&ENVIRON\n')
                for k, v in sorted(six.iteritems(environ_namelist)):
                    handle.write(convert_input_to_namelist_entry(k, v, mapping=mapping_species))
                handle.write('/\n')

        # Check for the deprecated 'ALSO_BANDS' setting and if present fire a deprecation log message
        also_bands = settings.pop('ALSO_BANDS', None)
        if also_bands:
            import logging
            from aiida.common.log import get_dblogger_extra

            logger = logging.LoggerAdapter(logger=self.logger, extra=get_dblogger_extra(self))
            logger.warning(
                "The '{}' setting is deprecated as bands are now parsed by default. "
                "If you do not want the bands to be parsed set the '{}' to True {}. "
                "Note that the eigenvalue.xml files are also no longer stored in the repository"
                .format('also_bands', 'no_bands', type(self))
            )

        calcinfo = datastructures.CalcInfo()

        calcinfo.uuid = str(self.uuid)
        # Empty command line by default
        cmdline_params = settings.pop('CMDLINE', [])
        # we commented calcinfo.stin_name and added it here in cmdline_params
        # in this way the mpirun ... pw.x ... < aiida.in
        # is replaced by mpirun ... pw.x ... -in aiida.in
        # in the scheduler, _get_run_line, if cmdline_params is empty, it
        # simply uses < calcinfo.stin_name
        calcinfo.cmdline_params = (list(cmdline_params) + ["-in", self._INPUT_FILE_NAME])

        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = (list(cmdline_params) + ["-in", self._INPUT_FILE_NAME])
        codeinfo.stdout_name = self._OUTPUT_FILE_NAME
        codeinfo.code_uuid = self.inputs.code.uuid
        calcinfo.codes_info = [codeinfo]

        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.remote_symlink_list = remote_symlink_list

        # Retrieve by default the output file and the xml file
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(self._OUTPUT_FILE_NAME)
        calcinfo.retrieve_list.extend(self.xml_filepaths)
        calcinfo.retrieve_list += settings.pop('ADDITIONAL_RETRIEVE_LIST', [])
        calcinfo.retrieve_list += self._internal_retrieve_list

        # Retrieve the k-point directories with the xml files to the temporary folder
        # to parse the band eigenvalues and occupations but not to have to save the raw files
        # if and only if the 'no_bands' key was not set to true in the settings
        no_bands = settings.pop('NO_BANDS', False)
        if no_bands is False:
            xmlpaths = os.path.join(self._OUTPUT_SUBFOLDER, self._PREFIX + '.save', 'K*[0-9]', 'eigenval*.xml')
            calcinfo.retrieve_temporary_list = [[xmlpaths, '.', 2]]

        # TODO: self.get_parserclass() probably raises AttributeError, but is caught later!
        #       Use self.input('metadata.options.parser_name') instead
        try:
            Parserclass = self.get_parserclass()
            parser = Parserclass(self)
            parser_opts = parser.get_parser_settings_key().upper()
            settings.pop(parser_opts)
        except (KeyError, AttributeError):
            # the key parser_opts isn't inside the dictionary
            pass

        if settings:
            unknown_keys = ', '.join(list(settings.keys()))
            raise exceptions.InputValidationError('`settings` contained unknown keys: {}'.format(unknown_keys))

        return calcinfo

    def _generate_PWCPinputdata(self, parameters, settings, pseudos, structure, kpoints=None):
        """
        This method creates the content of an input file
        in the PW/CP format.
        :
        """
        from aiida.common.utils import get_unique_filename
        import re
        local_copy_list_to_append = []

        # I put the first-level keys as uppercase (i.e., namelist and card names)
        # and the second-level keys as lowercase
        # (deeper levels are unchanged)
        input_params = _uppercase_dict(parameters.get_dict(), dict_name='parameters')
        input_params = {k: _lowercase_dict(v, dict_name=k) for k, v in six.iteritems(input_params)}

        # I remove unwanted elements (for the moment, instead, I stop; to change when
        # we setup a reasonable logging)
        for blocked in self._blocked_keywords:
            nl = blocked[0].upper()
            flag = blocked[1].lower()
            defaultvalue = None
            if len(blocked) >= 3:
                defaultvalue = blocked[2]
            if nl in input_params:
                # The following lines is meant to avoid putting in input the
                # parameters like celldm(*)
                stripped_inparams = [re.sub("[(0-9)]", "", _)
                                     for _ in input_params[nl].keys()]
                if flag in stripped_inparams:
                    raise exceptions.InputValidationError(
                        "You cannot specify explicitly the '{}' flag in the '{}' "
                        "namelist or card.".format(flag, nl))
                if defaultvalue is not None:
                    if nl not in input_params:
                        input_params[nl] = {}
                    input_params[nl][flag] = defaultvalue

        # Set some variables (look out at the case! NAMELISTS should be uppercase,
        # internal flag names must be lowercase)
        input_params.setdefault('CONTROL', {})
        input_params['CONTROL']['pseudo_dir'] = self._PSEUDO_SUBFOLDER
        input_params['CONTROL']['outdir'] = self._OUTPUT_SUBFOLDER
        input_params['CONTROL']['prefix'] = self._PREFIX
        input_params['CONTROL']['verbosity'] = input_params['CONTROL'].get('verbosity', self._default_verbosity)  # Set to high if not specified

        # ============ I prepare the input site data =============
        # ------------ CELL_PARAMETERS -----------
        cell_parameters_card = "CELL_PARAMETERS angstrom\n"
        for vector in structure.cell:
            cell_parameters_card += ("{0:18.10f} {1:18.10f} {2:18.10f}"
                                     "\n".format(*vector))

        # ------------- ATOMIC_SPECIES ------------
        atomic_species_card_list = []

        # Keep track of the filenames to avoid to overwrite files
        # I use a dictionary where the key is the pseudo PK and the value
        # is the filename I used. In this way, I also use the same filename
        # if more than one kind uses the same pseudo.
        pseudo_filenames = {}

        # I keep track of the order of species
        kind_names = []
        # I add the pseudopotential files to the list of files to be copied
        for kind in structure.kinds:
            # This should not give errors, I already checked before that
            # the list of keys of pseudos and kinds coincides
            ps = pseudos[kind.name]
            if kind.is_alloy or kind.has_vacancies:
                raise exceptions.InputValidationError("Kind '{}' is an alloy or has "
                                           "vacancies. This is not allowed for pw.x input structures."
                                           "".format(kind.name))

            try:
                # It it is the same pseudopotential file, use the same filename
                filename = pseudo_filenames[ps.pk]
            except KeyError:
                # The pseudo was not encountered yet; use a new name and also add it to the local copy list
                filename = get_unique_filename(ps.filename, list(pseudo_filenames.values()))
                pseudo_filenames[ps.pk] = filename
                local_copy_list_to_append.append((ps.uuid, ps.filename, os.path.join(self._PSEUDO_SUBFOLDER, filename)))

            kind_names.append(kind.name)
            atomic_species_card_list.append("{} {} {}\n".format(kind.name.ljust(6), kind.mass, filename))

        # I join the lines, but I resort them using the alphabetical order of
        # species, given by the kind_names list. I also store the mapping_species
        # list, with the order of species used in the file
        mapping_species, sorted_atomic_species_card_list = list(zip(
            *sorted(zip(kind_names, atomic_species_card_list))))
        # The format of mapping_species required later is a dictionary, whose
        # values are the indices, so I convert to this format
        # Note the (idx+1) to convert to fortran 1-based lists
        mapping_species = {sp_name: (idx + 1) for idx, sp_name
                           in enumerate(mapping_species)}
        # I add the first line
        sorted_atomic_species_card_list = (["ATOMIC_SPECIES\n"] +
                                           list(
                                               sorted_atomic_species_card_list))
        atomic_species_card = "".join(sorted_atomic_species_card_list)
        # Free memory
        del sorted_atomic_species_card_list
        del atomic_species_card_list

        # ------------ ATOMIC_POSITIONS -----------
        atomic_positions_card_list = ["ATOMIC_POSITIONS angstrom\n"]

        # Check on validity of FIXED_COORDS
        fixed_coords_strings = []
        fixed_coords = settings.pop('FIXED_COORDS', None)
        if fixed_coords is None:
            # No fixed_coords specified: I store a list of empty strings
            fixed_coords_strings = [""] * len(structure.sites)
        else:
            if len(fixed_coords) != len(structure.sites):
                raise exceptions.InputValidationError(
                    "Input structure contains {:d} sites, but "
                    "fixed_coords has length {:d}".format(len(structure.sites),
                                                          len(fixed_coords)))

            for i, this_atom_fix in enumerate(fixed_coords):
                if len(this_atom_fix) != 3:
                    raise exceptions.InputValidationError(
                        "fixed_coords({:d}) has not length three"
                        "".format(i + 1))
                for fixed_c in this_atom_fix:
                    if not isinstance(fixed_c, bool):
                        raise exceptions.InputValidationError(
                            "fixed_coords({:d}) has non-boolean "
                            "elements".format(i + 1))

                if_pos_values = [self._if_pos(_) for _ in this_atom_fix]
                fixed_coords_strings.append(
                    "  {:d} {:d} {:d}".format(*if_pos_values))

        for site, fixed_coords_string in zip(
                structure.sites, fixed_coords_strings):
            atomic_positions_card_list.append(
                "{0} {1:18.10f} {2:18.10f} {3:18.10f} {4}\n".format(
                    site.kind_name.ljust(6), site.position[0], site.position[1],
                    site.position[2], fixed_coords_string))
        atomic_positions_card = "".join(atomic_positions_card_list)
        del atomic_positions_card_list

        # Optional ATOMIC_FORCES card
        atomic_forces = settings.pop('ATOMIC_FORCES', None)
        if atomic_forces is not None:

            # Checking that there are as many forces defined as there are sites in the structure
            if len(atomic_forces) != len(structure.sites):
                raise exceptions.InputValidationError(
                    'Input structure contains {:d} sites, but atomic forces has length {:d}'.format(
                        len(structure.sites), len(atomic_forces)
                    )
                )

            lines = ['ATOMIC_FORCES\n']
            for site, vector in zip(structure.sites, atomic_forces):

                # Checking that all 3 dimensions are specified:
                if len(vector) != 3:
                    raise exceptions.InputValidationError('Forces({}) for {} has not length three'.format(vector, site))

                lines.append('{0} {1:18.10f} {2:18.10f} {3:18.10f}\n'.format(site.kind_name.ljust(6), *vector))

            # Append to atomic_positions_card so that this card will be printed directly after
            atomic_positions_card += ''.join(lines)
            del lines

        # Optional ATOMIC_VELOCITIES card
        atomic_velocities = settings.pop('ATOMIC_VELOCITIES', None)
        if atomic_velocities is not None:

            # Checking that there are as many velocities defined as there are sites in the structure
            if len(atomic_velocities) != len(structure.sites):
                raise exceptions.InputValidationError(
                    'Input structure contains {:d} sites, but atomic velocities has length {:d}'.format(
                        len(structure.sites), len(atomic_velocities)
                    )
                )

            lines = ['ATOMIC_VELOCITIES\n']
            for site, vector in zip(structure.sites, atomic_velocities):

                # Checking that all 3 dimensions are specified:
                if len(vector) != 3:
                    raise exceptions.InputValidationError('Velocities({}) for {} has not length three'.format(vector, site))

                lines.append('{0} {1:18.10f} {2:18.10f} {3:18.10f}\n'.format(site.kind_name.ljust(6), *vector))

            # Append to atomic_positions_card so that this card will be printed directly after
            atomic_positions_card += ''.join(lines)
            del lines

        # I set the variables that must be specified, related to the system
        # Set some variables (look out at the case! NAMELISTS should be
        # uppercase, internal flag names must be lowercase)
        input_params.setdefault('SYSTEM', {})
        input_params['SYSTEM']['ibrav'] = 0
        input_params['SYSTEM']['nat'] = len(structure.sites)
        input_params['SYSTEM']['ntyp'] = len(structure.kinds)

        # ============ I prepare the k-points =============
        if self._use_kpoints:
            try:
                mesh, offset = kpoints.get_kpoints_mesh()
                has_mesh = True
                force_kpoints_list = settings.pop('FORCE_KPOINTS_LIST', False)
                if force_kpoints_list:
                    kpoints_list = kpoints.get_kpoints_mesh(print_list=True)
                    num_kpoints = len(kpoints_list)
                    has_mesh = False
                    weights = [1.] * num_kpoints

            except AttributeError:

                try:
                    kpoints_list = kpoints.get_kpoints()
                    num_kpoints = len(kpoints_list)
                    has_mesh = False
                    if num_kpoints == 0:
                        raise exceptions.InputValidationError(
                            "At least one k point must be "
                            "provided for non-gamma calculations")
                except AttributeError:
                    raise exceptions.InputValidationError(
                        "No valid kpoints have been found")

                try:
                    _, weights = kpoints.get_kpoints(also_weights=True)
                except AttributeError:
                    weights = [1.] * num_kpoints

            gamma_only = settings.pop("GAMMA_ONLY", False)

            if gamma_only:
                if has_mesh:
                    if tuple(mesh) != (1, 1, 1) or tuple(offset) != (
                            0., 0., 0.):
                        raise exceptions.InputValidationError(
                            "If a gamma_only calculation is requested, the "
                            "kpoint mesh must be (1,1,1),offset=(0.,0.,0.)")

                else:
                    if (len(kpoints_list) != 1 or tuple(kpoints_list[0]) != tuple(0., 0., 0.)):
                        raise exceptions.InputValidationError(
                            "If a gamma_only calculation is requested, the "
                            "kpoints coordinates must only be (0.,0.,0.)")

                kpoints_type = "gamma"

            elif has_mesh:
                kpoints_type = "automatic"

            else:
                kpoints_type = "crystal"

            kpoints_card_list = ["K_POINTS {}\n".format(kpoints_type)]

            if kpoints_type == "automatic":
                if any([(i != 0. and i != 0.5) for i in offset]):
                    raise exceptions.InputValidationError("offset list must only be made "
                                               "of 0 or 0.5 floats")
                the_offset = [0 if i == 0. else 1 for i in offset]
                the_6_integers = list(mesh) + the_offset
                kpoints_card_list.append("{:d} {:d} {:d} {:d} {:d} {:d}\n"
                                         "".format(*the_6_integers))

            elif kpoints_type == "gamma":
                # nothing to be written in this case
                pass
            else:
                kpoints_card_list.append("{:d}\n".format(num_kpoints))
                for kpoint, weight in zip(kpoints_list, weights):
                    kpoints_card_list.append(
                        "  {:18.10f} {:18.10f} {:18.10f} {:18.10f}"
                        "\n".format(kpoint[0], kpoint[1], kpoint[2], weight))

            kpoints_card = "".join(kpoints_card_list)
            del kpoints_card_list

        # =================== NAMELISTS AND CARDS ========================
        try:
            namelists_toprint = settings.pop('NAMELISTS')
            if not isinstance(namelists_toprint, list):
                raise exceptions.InputValidationError(
                    "The 'NAMELISTS' value, if specified in the settings input "
                    "node, must be a list of strings")
        except KeyError:  # list of namelists not specified; do automatic detection
            try:
                control_nl = input_params['CONTROL']
                calculation_type = control_nl['calculation']
            except KeyError:
                raise exceptions.InputValidationError(
                    "No 'calculation' in CONTROL namelist."
                    "It is required for automatic detection of the valid list "
                    "of namelists. Otherwise, specify the list of namelists "
                    "using the NAMELISTS key inside the 'settings' input node")

            try:
                namelists_toprint = self._automatic_namelists[calculation_type]
            except KeyError:
                raise exceptions.InputValidationError("Unknown 'calculation' value in "
                                           "CONTROL namelist {}. Otherwise, specify the list of "
                                           "namelists using the NAMELISTS inside the 'settings' input "
                                           "node".format(calculation_type))

        inputfile = u''
        for namelist_name in namelists_toprint:
            inputfile += u'&{0}\n'.format(namelist_name)
            # namelist content; set to {} if not present, so that we leave an empty namelist
            namelist = input_params.pop(namelist_name, {})
            for key, value in sorted(namelist.items()):
                inputfile += convert_input_to_namelist_entry(key, value, mapping=mapping_species)
            inputfile += u'/\n'

        # Write cards now
        inputfile += atomic_species_card
        inputfile += atomic_positions_card
        if self._use_kpoints:
            inputfile += kpoints_card
        inputfile += cell_parameters_card

        if input_params:
            raise exceptions.InputValidationError(
                "The following namelists are specified in input_params, but are "
                "not valid namelists for the current type of calculation: "
                "{}".format(",".join(list(input_params.keys()))))

        return inputfile, local_copy_list_to_append

    def _if_pos(self, fixed):
        """
        Simple function that returns 0 if fixed is True, 1 otherwise.
        Useful to convert from the boolean value of fixed_coords to the value required
        by Quantum Espresso as if_pos.
        """
        if fixed:
            return 0
        else:
            return 1


def _lowercase_dict(d, dict_name):
    return _case_transform_dict(d, dict_name, "_lowercase_dict", str.lower)


def _uppercase_dict(d, dict_name):
    return _case_transform_dict(d, dict_name, "_uppercase_dict", str.upper)


def _case_transform_dict(d, dict_name, func_name, transform):
    from collections import Counter

    if not isinstance(d, dict):
        raise TypeError("{} accepts only dictionaries as argument, while you gave {}".format(func_name, type(d)))
    new_dict = dict((transform(str(k)), v) for k, v in six.iteritems(d))
    if len(new_dict) != len(d):
        num_items = Counter(transform(str(k)) for k in d.keys())
        double_keys = ",".join([k for k, v in num_items if v > 1])
        raise exceptions.InputValidationError(
            "Inside the dictionary '{}' there are the following keys that "
            "are repeated more than once when compared case-insensitively: {}."
            "This is not allowed.".format(dict_name, double_keys))
    return new_dict
