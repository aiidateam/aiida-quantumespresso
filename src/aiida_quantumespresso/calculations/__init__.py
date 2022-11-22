# -*- coding: utf-8 -*-
"""Base `CalcJob` for implementations for pw.x and cp.x of Quantum ESPRESSO."""
import abc
import copy
from functools import partial
import numbers
import os
from types import MappingProxyType
import warnings

from aiida import orm
from aiida.common import datastructures, exceptions
from aiida.common.lang import classproperty
from aiida.common.warnings import AiidaDeprecationWarning
from aiida.plugins import DataFactory
from qe_tools.converters import get_parameters_from_cell

from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry

from .base import CalcJob
from .helpers import QEInputValidationError

LegacyUpfData = DataFactory('core.upf')
UpfData = DataFactory('pseudo.upf')


class BasePwCpInputGenerator(CalcJob):
    """Base `CalcJob` for implementations for pw.x and cp.x of Quantum ESPRESSO."""

    _PSEUDO_SUBFOLDER = './pseudo/'
    _OUTPUT_SUBFOLDER = './out/'
    _PREFIX = 'aiida'
    _DEFAULT_INPUT_FILE = 'aiida.in'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'
    _DATAFILE_XML_PRE_6_2 = 'data-file.xml'
    _DATAFILE_XML_POST_6_2 = 'data-file-schema.xml'
    _ENVIRON_INPUT_FILE_NAME = 'environ.in'
    _DEFAULT_IBRAV = 0

    # A mapping {flag_name: help_string} of parallelization flags
    # possible in QE codes. The flags that are actually implemented in a
    # given code should be specified in the '_ENABLED_PARALLELIZATION_FLAGS'
    # tuple of each calculation subclass.
    _PARALLELIZATION_FLAGS = MappingProxyType(
        dict(
            nimage="The number of 'images', each corresponding to a different self-consistent or "
            'linear-response calculation.',
            npool="The number of 'pools', each taking care of a group of k-points.",
            nband="The number of 'band groups', each taking care of a group of Kohn-Sham orbitals.",
            ntg="The number of 'task groups' across which the FFT planes are distributed.",
            ndiag="The number of 'linear algebra groups' used when parallelizing the subspace "
            'diagonalization / iterative orthonormalization. By default, no parameter is '
            'passed to Quantum ESPRESSO, meaning it will use its default.',
            nhw="The 'nmany' FFT bands parallelization option."
        )
    )

    _ENABLED_PARALLELIZATION_FLAGS = tuple()

    _PARALLELIZATION_FLAG_ALIASES = MappingProxyType(
        dict(
            nimage=('ni', 'nimages', 'npot'),
            npool=('nk', 'npools'),
            nband=('nb', 'nbgrp', 'nband_group'),
            ntg=('nt', 'ntask_groups', 'nyfft'),
            ndiag=('northo', 'nd', 'nproc_diag', 'nproc_ortho'),
            nhw=('nh', 'n_howmany', 'howmany')
        )
    )

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
        # pylint: disable=no-self-argument
        return [cls._DATAFILE_XML_POST_6_2, cls._DATAFILE_XML_PRE_6_2]

    @abc.abstractmethod
    @classproperty
    def xml_filepaths(cls):  # pylint: disable=no-self-argument
        """Return a list of XML output filepaths relative to the remote working directory that should be retrieved."""

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('metadata.options.input_filename', valid_type=str, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=str, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.withmpi', valid_type=bool, default=True)  # Override default withmpi=False
        spec.input('structure', valid_type=orm.StructureData,
            help='The input structure.')
        spec.input('parameters', valid_type=orm.Dict,
            help='The input parameters that are to be used to construct the input file.')
        spec.input('settings', valid_type=orm.Dict, required=False,
            help='Optional parameters to affect the way the calculation job and the parsing are performed.')
        spec.input('parent_folder', valid_type=orm.RemoteData, required=False,
            help='An optional working directory of a previously completed calculation to restart from.')
        spec.input('vdw_table', valid_type=orm.SinglefileData, required=False,
            help='Optional van der Waals table contained in a `SinglefileData`.')
        spec.input_namespace('pseudos', valid_type=(LegacyUpfData, UpfData), dynamic=True, required=True,
            help='A mapping of `UpfData` nodes onto the kind name to which they should apply.')
        # yapf: enable
        spec.input(
            'parallelization',
            valid_type=orm.Dict,
            required=False,
            help=(
                'Parallelization options. The following flags are allowed:\n' + '\n'.join(
                    f'{flag_name:<7}: {cls._PARALLELIZATION_FLAGS[flag_name]}'
                    for flag_name in cls._ENABLED_PARALLELIZATION_FLAGS
                )
            ),
            validator=cls.validate_parallelization
        )
        spec.inputs.validator = cls.validate_inputs

    @classmethod
    def validate_inputs(cls, value, port_namespace):
        """Validate the entire inputs namespace."""

        # Wrapping processes may choose to exclude certain input ports in which case we can't validate. If the ports
        # have been excluded, and so are no longer part of the ``port_namespace``, skip the validation.
        if any(key not in port_namespace for key in ('pseudos', 'structure')):
            return

        # At this point, both ports are part of the namespace, and both are required so return an error message if any
        # of the two is missing.
        for key in ('pseudos', 'structure'):
            if key not in value:
                return f'required value was not provided for the `{key}` namespace.'

        structure_kinds = set(value['structure'].get_kind_names())
        pseudo_kinds = set(value['pseudos'].keys())

        if structure_kinds != pseudo_kinds:
            return f'The `pseudos` specified and structure kinds do not match: {pseudo_kinds} vs {structure_kinds}'

    @classmethod
    def validate_parallelization(cls, value, _):
        """Validate the ``parallelization`` port."""
        if value:
            value_dict = value.get_dict()
            unknown_flags = set(value_dict.keys()) - set(cls._ENABLED_PARALLELIZATION_FLAGS)
            if unknown_flags:
                return (
                    f"Unknown flags in 'parallelization': {unknown_flags}, "
                    f'allowed flags are {cls._ENABLED_PARALLELIZATION_FLAGS}.'
                )
            invalid_values = [val for val in value_dict.values() if not isinstance(val, numbers.Integral)]
            if invalid_values:
                return f'Parallelization values must be integers; got invalid values {invalid_values}.'

    def prepare_for_submission(self, folder):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        # pylint: disable=too-many-branches,too-many-statements
        if 'settings' in self.inputs:
            settings = _uppercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
        else:
            settings = {}

        # Check that a pseudo potential was specified for each kind present in the `StructureData`
        kinds = [kind.name for kind in self.inputs.structure.kinds]
        if set(kinds) != set(self.inputs.pseudos.keys()):
            formatted_pseudos = ', '.join(list(self.inputs.pseudos.keys()))
            formatted_kinds = ', '.join(list(kinds))
            raise exceptions.InputValidationError(
                'Mismatch between the defined pseudos and the list of kinds of the structure.\n'
                f'Pseudos: {formatted_pseudos};\nKinds: {formatted_kinds}'
            )

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
            uuid = self.inputs.vdw_table.uuid
            src_path = self.inputs.vdw_table.filename
            dst_path = os.path.join(self._PSEUDO_SUBFOLDER, self.inputs.vdw_table.filename)
            local_copy_list.append((uuid, src_path, dst_path))

        if 'hubbard_file' in self.inputs:
            uuid = self.inputs.hubbard_file.uuid
            src_path = self.inputs.hubbard_file.filename
            dst_path = self.filename_input_hubbard_parameters
            local_copy_list.append((uuid, src_path, dst_path))

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

        with folder.open(self.metadata.options.input_filename, 'w') as handle:
            handle.write(input_filecontent)

        # operations for restart
        symlink = settings.pop('PARENT_FOLDER_SYMLINK', self._default_symlink_usage)  # a boolean
        if symlink:
            if 'parent_folder' in self.inputs:
                # I put the symlink to the old parent ./out folder
                remote_symlink_list.append((
                    self.inputs.parent_folder.computer.uuid,
                    os.path.join(self.inputs.parent_folder.get_remote_path(),
                                 self._restart_copy_from), self._restart_copy_to
                ))
        else:
            # copy remote output dir, if specified
            if 'parent_folder' in self.inputs:
                remote_copy_list.append((
                    self.inputs.parent_folder.computer.uuid,
                    os.path.join(self.inputs.parent_folder.get_remote_path(),
                                 self._restart_copy_from), self._restart_copy_to
                ))

        # Create an `.EXIT` file if `only_initialization` flag in `settings` is set to `True`
        if settings.pop('ONLY_INITIALIZATION', False):
            with folder.open(f'{self._PREFIX}.EXIT', 'w') as handle:
                handle.write('\n')

        # Check if specific inputs for the ENVIRON module where specified
        environ_namelist = settings.pop('ENVIRON', None)
        if environ_namelist is not None:
            if not isinstance(environ_namelist, dict):
                raise exceptions.InputValidationError('ENVIRON namelist should be specified as a dictionary')
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

            with folder.open(self._ENVIRON_INPUT_FILE_NAME, 'w') as handle:
                handle.write('&ENVIRON\n')
                for key, value in sorted(environ_namelist.items()):
                    handle.write(convert_input_to_namelist_entry(key, value, mapping=mapping_species))
                handle.write('/\n')

        # Check for the deprecated 'ALSO_BANDS' setting and if present fire a deprecation log message
        also_bands = settings.pop('ALSO_BANDS', None)
        if also_bands:
            self.node.logger.warning(
                'The `also_bands` setting is deprecated as bands are now parsed by default. '
                'If you do not want the bands to be parsed set the `no_bands` to True. '
                'Note that the eigenvalue.xml files are also no longer stored in the repository'
            )

        calcinfo = datastructures.CalcInfo()

        calcinfo.uuid = str(self.uuid)
        # Start from an empty command line by default
        cmdline_params = self._add_parallelization_flags_to_cmdline_params(cmdline_params=settings.pop('CMDLINE', []))

        # we commented calcinfo.stin_name and added it here in cmdline_params
        # in this way the mpirun ... pw.x ... < aiida.in
        # is replaced by mpirun ... pw.x ... -in aiida.in
        # in the scheduler, _get_run_line, if cmdline_params is empty, it
        # simply uses < calcinfo.stin_name
        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = (list(cmdline_params) + ['-in', self.metadata.options.input_filename])
        codeinfo.stdout_name = self.metadata.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid
        calcinfo.codes_info = [codeinfo]

        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.remote_symlink_list = remote_symlink_list

        # Retrieve by default the output file and the xml file
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(self.metadata.options.output_filename)
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

        # We might still have parser options in the settings dictionary: pop them.
        _pop_parser_options(self, settings)

        if settings:
            unknown_keys = ', '.join(list(settings.keys()))
            raise exceptions.InputValidationError(f'`settings` contained unexpected keys: {unknown_keys}')

        return calcinfo

    def _add_parallelization_flags_to_cmdline_params(self, cmdline_params):
        """Get the command line parameters with added parallelization flags.

        Adds the parallelization flags to the given `cmdline_params` and
        returns the updated list.

        Raises an `InputValidationError` if multiple aliases to the same
        flag are given in `cmdline_params`, or the same flag is given
        both in `cmdline_params` and the explicit `parallelization`
        input.
        """
        cmdline_params_res = copy.deepcopy(cmdline_params)
        # The `cmdline_params_normalized` are used only here to check
        # for existing parallelization flags.
        cmdline_params_normalized = []
        for param in cmdline_params:
            cmdline_params_normalized.extend(param.split())

        if 'parallelization' in self.inputs:
            parallelization_dict = self.inputs.parallelization.get_dict()
        else:
            parallelization_dict = {}
        # To make the order of flags consistent and "nice", we use the
        # ordering from the flag definition.
        for flag_name in self._ENABLED_PARALLELIZATION_FLAGS:
            all_aliases = list(self._PARALLELIZATION_FLAG_ALIASES[flag_name]) + [flag_name]
            aliases_in_cmdline = [alias for alias in all_aliases if f'-{alias}' in cmdline_params_normalized]
            if aliases_in_cmdline:
                if len(aliases_in_cmdline) > 1:
                    raise exceptions.InputValidationError(
                        f'Conflicting parallelization flags {aliases_in_cmdline} '
                        "in settings['CMDLINE']"
                    )
                if flag_name in parallelization_dict:
                    raise exceptions.InputValidationError(
                        f"Parallelization flag '{aliases_in_cmdline[0]}' specified in settings['CMDLINE'] conflicts "
                        f"with '{flag_name}' in the 'parallelization' input."
                    )
                else:
                    warnings.warn(
                        "Specifying the parallelization flags through settings['CMDLINE'] is "
                        "deprecated, use the 'parallelization' input instead.", AiidaDeprecationWarning
                    )
                    continue
            if flag_name in parallelization_dict:
                flag_value = parallelization_dict[flag_name]
                cmdline_params_res += [f'-{flag_name}', str(flag_value)]
        return cmdline_params_res

    @staticmethod
    def _generate_PWCP_input_tail(*args, **kwargs):
        """Generate tail of input file.

        By default, nothing specific is generated.
        This method can be implemented again in derived classes, and it will be called by _generate_PWCPinputdata
        """
        # pylint: disable=unused-argument,invalid-name
        return ''

    @classmethod
    def _generate_PWCPinputdata(cls, parameters, settings, pseudos, structure, kpoints=None, use_fractional=False):  # pylint: disable=invalid-name
        """Create the input file in string format for a pw.x or cp.x calculation for the given inputs."""
        # pylint: disable=too-many-branches,too-many-statements
        import re

        from aiida.common.utils import get_unique_filename
        local_copy_list_to_append = []

        # I put the first-level keys as uppercase (i.e., namelist and card names)
        # and the second-level keys as lowercase
        # (deeper levels are unchanged)
        input_params = _uppercase_dict(parameters.get_dict(), dict_name='parameters')
        input_params = {k: _lowercase_dict(v, dict_name=k) for k, v in input_params.items()}

        # I remove unwanted elements (for the moment, instead, I stop; to change when we setup a reasonable logging)
        for blocked in cls._blocked_keywords:
            namelist = blocked[0].upper()
            flag = blocked[1].lower()
            defaultvalue = None
            if len(blocked) >= 3:
                defaultvalue = blocked[2]
            if namelist in input_params:
                # The following lines is meant to avoid putting in input the
                # parameters like celldm(*)
                stripped_inparams = [re.sub('[(0-9)]', '', _) for _ in input_params[namelist].keys()]
                if flag in stripped_inparams:
                    raise exceptions.InputValidationError(
                        f"You cannot specify explicitly the '{flag}' flag in the '{namelist}' namelist or card."
                    )
                if defaultvalue is not None:
                    if namelist not in input_params:
                        input_params[namelist] = {}
                    input_params[namelist][flag] = defaultvalue

        # Set some variables (look out at the case! NAMELISTS should be uppercase,
        # internal flag names must be lowercase)
        input_params.setdefault('CONTROL', {})
        input_params['CONTROL']['pseudo_dir'] = cls._PSEUDO_SUBFOLDER
        input_params['CONTROL']['outdir'] = cls._OUTPUT_SUBFOLDER
        input_params['CONTROL']['prefix'] = cls._PREFIX
        input_params['CONTROL']['verbosity'] = input_params['CONTROL'].get('verbosity', cls._default_verbosity)

        # ============ I prepare the input site data =============
        # ------------ CELL_PARAMETERS -----------

        # Specify cell parameters only if 'ibrav' is zero.
        if input_params.get('SYSTEM', {}).get('ibrav', cls._DEFAULT_IBRAV) == 0:
            cell_parameters_card = 'CELL_PARAMETERS angstrom\n'
            for vector in structure.cell:
                cell_parameters_card += ('{0:18.10f} {1:18.10f} {2:18.10f}\n'.format(*vector))  # pylint: disable=consider-using-f-string
        else:
            cell_parameters_card = ''

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
            pseudo = pseudos[kind.name]
            if kind.is_alloy or kind.has_vacancies:
                raise exceptions.InputValidationError(
                    f"Kind '{kind.name}' is an alloy or has vacancies. This is not allowed for pw.x input structures."
                )

            try:
                # If it is the same pseudopotential file, use the same filename
                filename = pseudo_filenames[pseudo.pk]
            except KeyError:
                # The pseudo was not encountered yet; use a new name and also add it to the local copy list
                filename = get_unique_filename(pseudo.filename, list(pseudo_filenames.values()))
                pseudo_filenames[pseudo.pk] = filename
                local_copy_list_to_append.append(
                    (pseudo.uuid, pseudo.filename, os.path.join(cls._PSEUDO_SUBFOLDER, filename))
                )

            kind_names.append(kind.name)
            atomic_species_card_list.append(f'{kind.name.ljust(6)} {kind.mass} {filename}\n')

        # I join the lines, but I resort them using the alphabetical order of
        # species, given by the kind_names list. I also store the mapping_species
        # list, with the order of species used in the file
        mapping_species, sorted_atomic_species_card_list = list(zip(*sorted(zip(kind_names, atomic_species_card_list))))
        # The format of mapping_species required later is a dictionary, whose
        # values are the indices, so I convert to this format
        # Note the (idx+1) to convert to fortran 1-based lists
        mapping_species = {sp_name: (idx + 1) for idx, sp_name in enumerate(mapping_species)}
        # I add the first line
        sorted_atomic_species_card_list = (['ATOMIC_SPECIES\n'] + list(sorted_atomic_species_card_list))
        atomic_species_card = ''.join(sorted_atomic_species_card_list)
        # Free memory
        del sorted_atomic_species_card_list
        del atomic_species_card_list

        # ------------ ATOMIC_POSITIONS -----------
        # Check on validity of FIXED_COORDS
        fixed_coords_strings = []
        fixed_coords = settings.pop('FIXED_COORDS', None)
        if fixed_coords is None:
            # No fixed_coords specified: I store a list of empty strings
            fixed_coords_strings = [''] * len(structure.sites)
        else:
            if len(fixed_coords) != len(structure.sites):
                raise exceptions.InputValidationError(
                    f'Input structure contains {len(structure.sites)} sites, but fixed_coords has length '
                    f'{len(fixed_coords)}'
                )

            for i, this_atom_fix in enumerate(fixed_coords):
                if len(this_atom_fix) != 3:
                    raise exceptions.InputValidationError(f'fixed_coords({i + 1:d}) has not length three')
                for fixed_c in this_atom_fix:
                    if not isinstance(fixed_c, bool):
                        raise exceptions.InputValidationError(f'fixed_coords({i + 1:d}) has non-boolean elements')

                if_pos_values = [cls._if_pos(_) for _ in this_atom_fix]
                fixed_coords_strings.append('  {:d} {:d} {:d}'.format(*if_pos_values))  # pylint: disable=consider-using-f-string

        abs_pos = [_.position for _ in structure.sites]
        if use_fractional:
            import numpy as np
            atomic_positions_card_list = ['ATOMIC_POSITIONS crystal\n']
            coordinates = np.dot(np.array(abs_pos), np.linalg.inv(np.array(structure.cell)))
        else:
            atomic_positions_card_list = ['ATOMIC_POSITIONS angstrom\n']
            coordinates = abs_pos

        for site, site_coords, fixed_coords_string in zip(structure.sites, coordinates, fixed_coords_strings):
            atomic_positions_card_list.append(
                '{0} {1:18.10f} {2:18.10f} {3:18.10f} {4}\n'.format(  # pylint: disable=consider-using-f-string
                    site.kind_name.ljust(6), site_coords[0], site_coords[1], site_coords[2], fixed_coords_string
                )
            )

        atomic_positions_card = ''.join(atomic_positions_card_list)
        del atomic_positions_card_list

        # Optional ATOMIC_FORCES card
        atomic_forces = settings.pop('ATOMIC_FORCES', None)
        if atomic_forces is not None:

            # Checking that there are as many forces defined as there are sites in the structure
            if len(atomic_forces) != len(structure.sites):
                raise exceptions.InputValidationError(
                    f'Input structure contains {len(structure.sites):d} sites, but atomic forces has length '
                    f'{len(atomic_forces):d}'
                )

            lines = ['ATOMIC_FORCES\n']
            for site, vector in zip(structure.sites, atomic_forces):

                # Checking that all 3 dimensions are specified:
                if len(vector) != 3:
                    raise exceptions.InputValidationError(f'Forces({vector}) for {site} has not length three')

                lines.append('{0} {1:18.10f} {2:18.10f} {3:18.10f}\n'.format(site.kind_name.ljust(6), *vector))  # pylint: disable=consider-using-f-string

            # Append to atomic_positions_card so that this card will be printed directly after
            atomic_positions_card += ''.join(lines)
            del lines

        # Optional ATOMIC_VELOCITIES card
        atomic_velocities = settings.pop('ATOMIC_VELOCITIES', None)
        if atomic_velocities is not None:

            # Checking that there are as many velocities defined as there are sites in the structure
            if len(atomic_velocities) != len(structure.sites):
                raise exceptions.InputValidationError(
                    f'Input structure contains {len(structure.sites):d} sites, but atomic velocities has length '
                    f'{len(atomic_velocities):d}'
                )

            lines = ['ATOMIC_VELOCITIES\n']
            for site, vector in zip(structure.sites, atomic_velocities):

                # Checking that all 3 dimensions are specified:
                if len(vector) != 3:
                    raise exceptions.InputValidationError(f'Velocities({vector}) for {site} has not length three')

                lines.append('{0} {1:18.10f} {2:18.10f} {3:18.10f}\n'.format(site.kind_name.ljust(6), *vector))  # pylint: disable=consider-using-f-string

            # Append to atomic_positions_card so that this card will be printed directly after
            atomic_positions_card += ''.join(lines)
            del lines

        # I set the variables that must be specified, related to the system
        # Set some variables (look out at the case! NAMELISTS should be
        # uppercase, internal flag names must be lowercase)
        input_params.setdefault('SYSTEM', {})
        input_params['SYSTEM'].setdefault('ibrav', cls._DEFAULT_IBRAV)
        ibrav = input_params['SYSTEM']['ibrav']
        if ibrav != 0:
            try:
                structure_parameters = get_parameters_from_cell(
                    ibrav=ibrav,
                    cell=structure.base.attributes.get('cell'),
                    tolerance=settings.pop('IBRAV_CELL_TOLERANCE', 1e-6)
                )
            except ValueError as exc:
                raise QEInputValidationError(f'Cannot get structure parameters from cell: {exc}') from exc
            input_params['SYSTEM'].update(structure_parameters)
        input_params['SYSTEM']['nat'] = len(structure.sites)
        input_params['SYSTEM']['ntyp'] = len(structure.kinds)

        # ============ I prepare the k-points =============
        if cls._use_kpoints:
            try:
                mesh, offset = kpoints.get_kpoints_mesh()
                has_mesh = True
                force_kpoints_list = settings.pop('FORCE_KPOINTS_LIST', False)
                if force_kpoints_list:
                    kpoints_list = kpoints.get_kpoints_mesh(print_list=True)
                    num_kpoints = len(kpoints_list)
                    has_mesh = False
                    weights = [1.] * num_kpoints

            except AttributeError as exception:

                try:
                    kpoints_list = kpoints.get_kpoints()
                    num_kpoints = len(kpoints_list)
                    has_mesh = False
                    if num_kpoints == 0:
                        raise exceptions.InputValidationError(
                            'At least one k point must be provided for non-gamma calculations'
                        ) from exception
                except AttributeError:
                    raise exceptions.InputValidationError('No valid kpoints have been found') from exception

                try:
                    _, weights = kpoints.get_kpoints(also_weights=True)
                except AttributeError:
                    weights = [1.] * num_kpoints

            gamma_only = settings.pop('GAMMA_ONLY', False)

            if gamma_only:
                if has_mesh:
                    if tuple(mesh) != (1, 1, 1) or tuple(offset) != (0., 0., 0.):
                        raise exceptions.InputValidationError(
                            'If a gamma_only calculation is requested, the '
                            'kpoint mesh must be (1,1,1),offset=(0.,0.,0.)'
                        )

                else:
                    if (len(kpoints_list) != 1 or tuple(kpoints_list[0]) != tuple(0., 0., 0.)):
                        raise exceptions.InputValidationError(
                            'If a gamma_only calculation is requested, the '
                            'kpoints coordinates must only be (0.,0.,0.)'
                        )

                kpoints_type = 'gamma'

            elif has_mesh:
                kpoints_type = 'automatic'

            else:
                kpoints_type = 'crystal'

            kpoints_card_list = [f'K_POINTS {kpoints_type}\n']

            if kpoints_type == 'automatic':
                if any(i not in [0, 0.5] for i in offset):
                    raise exceptions.InputValidationError('offset list must only be made of 0 or 0.5 floats')
                the_offset = [0 if i == 0. else 1 for i in offset]
                the_6_integers = list(mesh) + the_offset
                kpoints_card_list.append('{:d} {:d} {:d} {:d} {:d} {:d}\n'.format(*the_6_integers))  # pylint: disable=consider-using-f-string

            elif kpoints_type == 'gamma':
                # nothing to be written in this case
                pass
            else:
                kpoints_card_list.append(f'{num_kpoints:d}\n')
                for kpoint, weight in zip(kpoints_list, weights):
                    kpoints_card_list.append(
                        f'  {kpoint[0]:18.10f} {kpoint[1]:18.10f} {kpoint[2]:18.10f} {weight:18.10f}\n'
                    )

            kpoints_card = ''.join(kpoints_card_list)
            del kpoints_card_list

        # =================== NAMELISTS AND CARDS ========================
        try:
            namelists_toprint = settings.pop('NAMELISTS')
            if not isinstance(namelists_toprint, list):
                raise exceptions.InputValidationError(
                    "The 'NAMELISTS' value, if specified in the settings input "
                    'node, must be a list of strings'
                )
        except KeyError:  # list of namelists not specified; do automatic detection
            try:
                control_nl = input_params['CONTROL']
                calculation_type = control_nl['calculation']
            except KeyError as exception:
                raise exceptions.InputValidationError(
                    "No 'calculation' in CONTROL namelist."
                    'It is required for automatic detection of the valid list '
                    'of namelists. Otherwise, specify the list of namelists '
                    "using the NAMELISTS key inside the 'settings' input node."
                ) from exception

            try:
                namelists_toprint = cls._automatic_namelists[calculation_type]
            except KeyError as exception:
                raise exceptions.InputValidationError(
                    'Unknown `calculation` value in CONTROL namelist {calculation_type}. Otherwise, specify the list of'
                    'namelists using the NAMELISTS inside the `settings` input node'
                ) from exception

        inputfile = ''
        for namelist_name in namelists_toprint:
            inputfile += f'&{namelist_name}\n'
            # namelist content; set to {} if not present, so that we leave an empty namelist
            namelist = input_params.pop(namelist_name, {})
            for key, value in sorted(namelist.items()):
                inputfile += convert_input_to_namelist_entry(key, value, mapping=mapping_species)
            inputfile += '/\n'

        # Write cards now
        inputfile += atomic_species_card
        inputfile += atomic_positions_card
        if cls._use_kpoints:
            inputfile += kpoints_card
        inputfile += cell_parameters_card

        # Generate additional cards bases on input parameters and settings that are subclass specific
        tail = cls._generate_PWCP_input_tail(input_params=input_params, settings=settings)
        if tail:
            inputfile += f'\n{tail}'

        if input_params:
            raise exceptions.InputValidationError(
                'The following namelists are specified in input_params, but are not valid namelists for the current '
                f'type of calculation: {",".join(list(input_params.keys()))}'
            )

        return inputfile, local_copy_list_to_append

    @staticmethod
    def _if_pos(fixed):
        """Return 0 if fixed is True, 1 otherwise.

        Useful to convert from the boolean value of fixed_coords to the value required by Quantum Espresso as if_pos.
        """
        if fixed:
            return 0

        return 1


def _lowercase_dict(dictionary, dict_name):
    return _case_transform_dict(dictionary, dict_name, '_lowercase_dict', str.lower)


def _uppercase_dict(dictionary, dict_name):
    return _case_transform_dict(dictionary, dict_name, '_uppercase_dict', str.upper)


def _case_transform_dict(dictionary, dict_name, func_name, transform):
    from collections import Counter

    if not isinstance(dictionary, dict):
        raise TypeError(f'{func_name} accepts only dictionaries as argument, got {type(dictionary)}')
    new_dict = dict((transform(str(k)), v) for k, v in dictionary.items())
    if len(new_dict) != len(dictionary):
        num_items = Counter(transform(str(k)) for k in dictionary.keys())
        double_keys = ','.join([k for k, v in num_items if v > 1])
        raise exceptions.InputValidationError(
            f'Inside the dictionary `{dict_name}` there are the following keys that are repeated more than once when '
            f'compared case-insensitively: {double_keys}. This is not allowed.'
        )
    return new_dict


def _pop_parser_options(calc_job_instance, settings_dict, ignore_errors=True):
    """Delete any parser options from the settings dictionary.

    The parser options key is found via the get_parser_settings_key() method of the parser class specified as a metadata
    input.
    """
    from aiida.common import EntryPointError
    from aiida.plugins import ParserFactory
    try:
        parser_name = calc_job_instance.inputs['metadata']['options']['parser_name']
        parser_class = ParserFactory(parser_name)
        parser_opts_key = parser_class.get_parser_settings_key().upper()
        return settings_dict.pop(parser_opts_key, None)
    except (KeyError, EntryPointError, AttributeError) as exc:
        # KeyError: input 'metadata.options.parser_name' is not defined;
        # EntryPointError: there was an error loading the parser class form its entry point
        #   (this will probably cause errors elsewhere too);
        # AttributeError: the parser class doesn't have a method get_parser_settings_key().
        if ignore_errors:
            pass
        else:
            raise exc
