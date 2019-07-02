# -*- coding: utf-8 -*-
"""
Plugin to create a Quantum Espresso neb.x input file.
"""
from __future__ import absolute_import
import os
import copy
from aiida.common import InputValidationError
from aiida.common import CalcInfo, CodeInfo
from aiida.common import LinkType
from aiida.common.lang import classproperty
from aiida import orm
from aiida.engine import CalcJob
from aiida_quantumespresso.calculations import BasePwCpInputGenerator
from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.calculations import _lowercase_dict, _uppercase_dict
from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry
import six


class NebCalculation(BasePwCpInputGenerator, CalcJob):
    """
    Nudged Elastic Band code (neb.x) of Quantum ESPRESSO distribution
    For more information, refer to http://www.quantum-espresso.org/
    """

    _PREFIX = 'aiida'

    # in restarts, will not copy but use symlinks
    _default_symlink_usage = False

    # Default input and output file names
    # TODO: fix this mess: https://github.com/aiidateam/aiida-quantumespresso/issues/351
    _DEFAULT_INPUT_FILE = 'neb.dat'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'

    _automatic_namelists = {
        'scf': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
    }

    # Keywords that cannot be set (for the PW input)
    _blocked_keywords = [
        ('CONTROL', 'pseudo_dir'),  # set later
        ('CONTROL', 'outdir'),  # set later
        ('CONTROL', 'prefix'),  # set later
        ('SYSTEM', 'ibrav'),  # set later
        ('SYSTEM', 'celldm'),
        ('SYSTEM', 'nat'),  # set later
        ('SYSTEM', 'ntyp'),  # set later
        ('SYSTEM', 'a'), ('SYSTEM', 'b'), ('SYSTEM', 'c'),
        ('SYSTEM', 'cosab'), ('SYSTEM', 'cosac'), ('SYSTEM', 'cosbc'),
    ]

    # TODO: turn into input?
    _use_kpoints = True

    @classproperty
    def _internal_retrieve_list(cls):
        # I retrieve them all, even if I don't parse all of them
        _neb_ext_list = ['path', 'dat', 'int']
        return [ '{}.{}'.format(cls._PREFIX, ext) for ext in _neb_ext_list]

    @classproperty
    def xml_filepaths(cls):
        """Returns a list of relative filepaths of XML files."""
        filepaths = []

        for filename in cls.xml_filenames:
            filepath = os.path.join(cls._OUTPUT_SUBFOLDER, '{}_*[0-9].save'.format(cls._PREFIX), filename)
            filepaths.append(filepath)

        return filepaths

    @classmethod
    def define(cls, spec):
        CalcJob.define(spec)
        spec.input('metadata.options.input_filename', valid_type=six.string_types, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=six.string_types, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.parser_name', valid_type=six.string_types, default='quantumespresso.neb')
        spec.input('first_structure', valid_type=orm.StructureData, help='Initial structure')
        spec.input('last_structure', valid_type=orm.StructureData, help='Final structure')
        spec.input('parameters', valid_type=orm.Dict, help='NEB-specific input parameters')
        spec.input('settings', valid_type=orm.Dict, required=False,
            help='Optional parameters to affect the way the calculation job and the parsing are performed.')
        # We reuse some inputs from PwCalculation to construct the PW-specific parts of the input files:
        # 'parameters', 'pseudos', 'kpoints', ... TODO: what about settings?
        spec.expose_inputs(PwCalculation, namespace='pw', exclude=('structure','settings'))
        # TODO: outputs, exit codes

    # @classmethod
    # def define(cls, spec):
    #     super(NebCalculation, cls).define(spec)
    #     spec.input('metadata.options.input_filename', valid_type=six.string_types, default=cls._DEFAULT_INPUT_FILE)
    #     spec.input('metadata.options.output_filename', valid_type=six.string_types, default=cls._DEFAULT_OUTPUT_FILE)
    #     spec.input('metadata.options.parser_name', valid_type=six.string_types, default='quantumespresso.neb')
    #     spec.input('first_structure', valid_type=orm.StructureData, help='Choose the initial structure to use')  # TODO: better help strings
    #     spec.input('last_structure', valid_type=orm.StructureData, help='Choose the final structure to use')
    #     spec.input('kpoints', valid_type=orm.KpointsData, help='kpoint mesh or kpoint path to use')
    #     # TODO: how do I delete input 'parameters' from the base class?
    #     spec.input('pw_parameters', valid_type=orm.Dict,
    #         help='Input parameters used to construct the PW input file(s).')
    #     spec.input('neb_parameters', valid_type=orm.Dict,
    #         help='Input parameters used to construct the NEB input file.')
    #     # TODO: outputs, exit codes
    #     # spec.input('hubbard_file', valid_type=orm.SinglefileData, required=False,
    #     #     help='SinglefileData node containing the output Hubbard parameters from a HpCalculation')
    #     # spec.output('output_parameters', valid_type=orm.Dict,
    #     #     help='The `output_parameters` output node of the successful calculation.')
    #     # spec.output('output_structure', valid_type=orm.StructureData, required=False,
    #     #     help='The `output_structure` output node of the successful calculation if present.')
    #     # spec.output('output_trajectory', valid_type=orm.TrajectoryData, required=False)
    #     # spec.output('output_array', valid_type=orm.ArrayData, required=False,
    #     #     help='The `output_array` output node of the successful calculation if present.')
    #     # spec.output('output_band', valid_type=orm.BandsData, required=False,
    #     #     help='The `output_band` output node of the successful calculation if present.')
    #     # spec.output('output_kpoints', valid_type=orm.KpointsData, required=False)
    #     # spec.output('output_atomic_occupations', valid_type=orm.Dict, required=False)
    #     # spec.default_output_node = 'output_parameters'

    def _generate_NEBinputdata(self, neb_parameters, settings_dict):
        """ 
        This methods generate the input data for the NEB part of the calculation
        """
        # I put the first-level keys as uppercase (i.e., namelist and card names)
        # and the second-level keys as lowercase
        # (deeper levels are unchanged)
        input_params = _uppercase_dict(neb_parameters.get_dict(), dict_name='parameters')
        input_params = {k: _lowercase_dict(v, dict_name=k) for k, v in six.iteritems(input_params)}
        
        # For the neb input there are no blocked keywords
        
        # Create an empty dictionary for the compulsory namelist 'PATH' if not present
        if 'PATH' not in input_params:
            input_params['PATH'] = {}

        # In case of climbing image, we need the corresponding card
        climbing_image = False
        if input_params['PATH'].get('ci_scheme', 'no-ci').lower() in ['manual']:
            # TODO: why not 'auto' ?
            climbing_image = True
            try:
                climbing_image_list = settings_dict.pop("CLIMBING_IMAGES")
            except KeyError:
                raise InputValidationError("No climbing image specified for this calculation")
            if not isinstance(climbing_image_list, list):
                raise InputValidationError("Climbing images should be provided as a list")
            num_of_images = input_params['PATH'].get('num_of_images', 2)
            if any([ (i<2 or i>=num_of_images) for i in climbing_image_list ]):
                raise InputValidationError("The climbing images should be in the range between the first "
                                           "and the last image")

            climbing_image_card = "CLIMBING_IMAGES\n"
            climbing_image_card += ", ".join([str(_) for _ in climbing_image_list]) + "\n"
 
        inputfile = u"&PATH\n"
        # namelist content; set to {} if not present, so that we leave an empty namelist
        namelist = input_params.pop('PATH', {})
        for k, v in sorted(six.iteritems(namelist)):
            inputfile += convert_input_to_namelist_entry(k, v)
        inputfile += u"/\n"

        # Write cards now
        if climbing_image:
            inputfile += climbing_image_card

        if input_params:
            raise InputValidationError(
                "The following namelists are specified in input_params, but are "
                "not valid namelists for the current type of calculation: "
                "{}".format(",".join(list(input_params.keys()))))

        return inputfile

    def prepare_for_submission(self, folder):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        
        # TODO: call BasePwCpInputGenerator.prepare_for_submission (!?)
        
        import numpy as np

        local_copy_list = []
        remote_copy_list = []
        remote_symlink_list = []
        
        # Convert settings dictionary to have uppercase keys, or create an empty one if none was given.
        if 'settings' in self.inputs:
            settings_dict = _uppercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
        else:
            settings_dict = {}

        # TODO: remove these debug lines
        self.logger.debug('self.inputs.pw.pseudos has type <{}>', type(self.inputs.pw.pseudos))
        pseudos = self.inputs.pw.pseudos  # This is a PortNamespace, but can be used as a dict
        parent_calc_folder = self.inputs.get('parent_folder', None)
        vdw_table = self.inputs.get('vdw_table', None)

        # TODO: image != structure ?
        # Check that the first and last image have the same cell
        if abs(np.array(self.inputs.first_structure.cell)-
               np.array(self.inputs.last_structure.cell)).max() > 1.e-4:
            raise InputValidationError("Different cell in the fist and last image")

        # Check that the first and last image have the same number of sites
        if len(self.inputs.first_structure.sites) != len(self.inputs.last_structure.sites):
            raise InputValidationError("Different number of sites in the fist and last image")

        # Check that sites in the initial and final structure have the same kinds
        if self.inputs.first_structure.get_site_kindnames() != self.inputs.last_structure.get_site_kindnames():
            raise InputValidationError("Mismatch between the kind names and/or order between "
                                       "the first and final image")

        # Check that a pseudo potential was specified for each kind present in the `StructureData`
        kindnames = [kind.name for kind in self.inputs.first_structure.kinds]
        if set(kindnames) != set(self.inputs.pw.pseudos.keys()):
            raise InputValidationError(
                'Mismatch between the defined pseudos and the list of kinds of the structure.\n'
                'Pseudos: {};\nKinds: {}'.format(', '.join(list(self.inputs.pw.pseudos.keys())), ', '.join(list(kindnames))))

        ##############################
        # END OF INITIAL INPUT CHECK #
        ##############################

        # Create the subfolder that will contain the pseudopotentials
        folder.get_subfolder(self._PSEUDO_SUBFOLDER, create=True)
        # Create the subfolder for the output data (sometimes Quantum ESPRESSO codes crash if the folder does not exist)
        folder.get_subfolder(self._OUTPUT_SUBFOLDER, create=True)

        # We first prepare the NEB-specific input file.
        input_filecontent = self._generate_NEBinputdata(self.inputs.parameters, settings_dict)

        input_filename = folder.get_abs_path(self._INPUT_FILE_NAME)
        with open(input_filename, 'w') as handle:
            handle.write(input_filecontent)

        # We now generate the PW input files for each input structure
        local_copy_pseudo_list = []
        for i, structure in enumerate([self.inputs.first_structure, self.inputs.last_structure]):
            # We need to a pass a copy of the settings_dict for each structure
            this_settings_dict = copy.deepcopy(settings_dict)
            input_filecontent, this_local_copy_pseudo_list = self._generate_PWCPinputdata(
                self.inputs.pw.parameters, this_settings_dict, self.inputs.pw.pseudos, structure, self.inputs.pw.kpoints
            )
            local_copy_pseudo_list += this_local_copy_pseudo_list

            input_filename = folder.get_abs_path('pw_{}.in'.format(i+1))
            with open(input_filename, 'w') as handle:
                handle.write(input_filecontent)

        # We need to pop the settings that were used in the PW calculations
        for key in settings_dict.keys():
            if key not in list(this_settings_dict.keys()):
                settings_dict.pop(key)

        # We avoid to copy twice the same pseudopotential to the same filename
        local_copy_pseudo_list = set(local_copy_pseudo_list)
        # We check that two different pseudopotentials are not copied 
        # with the same name (otherwise the first is overwritten)
        if len({ filename for (uuid, filename, local_path) in local_copy_pseudo_list}) < len(local_copy_pseudo_list):
            raise InputValidationError("Same filename for two different pseudopotentials")

        local_copy_list += local_copy_pseudo_list 

        # If present, add also the Van der Waals table to the pseudo dir. Note that the name of the table is not checked
        # but should be the one expected by Quantum ESPRESSO.
        if vdw_table:
            local_copy_list.append((
                vdw_table.uuid,
                vdw_table.filename,
                os.path.join(self._PSEUDO_SUBFOLDER, vdw_table.filename)
            ))

        # operations for restart
        symlink = settings_dict.pop('PARENT_FOLDER_SYMLINK', self._default_symlink_usage)  # a boolean
        if symlink:
            if parent_calc_folder is not None:
                # I put the symlink to the old parent ./out folder
                remote_symlink_list.append((
                    parent_calc_folder.computer.uuid,
                    os.path.join(parent_calc_folder.get_remote_path(),
                                 self._OUTPUT_SUBFOLDER, '*'),  # asterisk: make individual symlinks for each file
                    self._OUTPUT_SUBFOLDER
                ))
                # and to the old parent prefix.path
                remote_symlink_list.append((
                    parent_calc_folder.computer.uuid,
                    os.path.join(parent_calc_folder.get_remote_path(),
                                 '{}.path'.format(self._PREFIX)),
                    '{}.path'.format(self._PREFIX)
                ))
        else:
            # copy remote output dir and .path file, if specified
            if parent_calc_folder is not None:
                remote_copy_list.append((
                    parent_calc_folder.computer.uuid,
                    os.path.join(parent_calc_folder.get_remote_path(),
                                 self._OUTPUT_SUBFOLDER, '*'),
                    self._OUTPUT_SUBFOLDER
                ))
                # and copy the old parent prefix.path
                remote_copy_list.append((
                    parent_calc_folder.computer.uuid,
                    os.path.join(parent_calc_folder.get_remote_path(),
                                 '{}.path'.format(self._PREFIX)),
                    '{}.path'.format(self._PREFIX)
                ))

        # here we may create an aiida.EXIT file
        create_exit_file = settings_dict.pop('ONLY_INITIALIZATION', False)
        if create_exit_file:
            exit_filename = '{}.EXIT'.format(self._PREFIX)
            with folder.open(exit_filename, 'w') as f:
                f.write('\n')

        calcinfo = CalcInfo()
        codeinfo = CodeInfo()

        calcinfo.uuid = self.uuid
        # Empty command line by default
        cmdline_params = settings_dict.pop('CMDLINE', [])
        # For the time-being we only have the initial and final image
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.remote_symlink_list = remote_symlink_list
        # In neb calculations there is no input read from standard input!! 

        codeinfo.cmdline_params = (["-input_images", "2"]
                                   + list(cmdline_params))
        codeinfo.stdout_name = self._OUTPUT_FILE_NAME
        codeinfo.code_uuid = self.inputs.code.uuid
        calcinfo.codes_info = [codeinfo]

        # Retrieve the output files and the xml files
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(self._OUTPUT_FILE_NAME)
        calcinfo.retrieve_list.append((
            os.path.join(self._OUTPUT_SUBFOLDER, self._PREFIX + '_*[0-9]', 'PW.out'),  # source relative path (globbing)
             '.',  # destination relative path
             2  # depth to preserve
        ))

        # TODO: check if anything has changed in the xml output (formats, names, contents)
        for xml_filepath in self.xml_filepaths:
            calcinfo.retrieve_list.append([xml_filepath, '.', 3])

        calcinfo.retrieve_list += settings_dict.pop('ADDITIONAL_RETRIEVE_LIST', [])
        calcinfo.retrieve_list += self._internal_retrieve_list

        # TODO: self.get_parserclass() probably raises AttributeError, but is caught later!
        #       Use self.input('metadata.options.parser_name') instead
        if settings_dict:
            try:
                Parserclass = self.get_parserclass()
                parser = Parserclass(self)
                parser_opts = parser.get_parser_settings_key()
                settings_dict.pop(parser_opts)
            except (KeyError, AttributeError):
                # the settings dictionary has no key 'parser_options',
                # or some methods don't exist
                pass

        if settings_dict:
            unknown_keys = ', '.join(list(settings_dict.keys()))
            raise InputValidationError('`settings` contained unknown keys: {}'.format(unknown_keys))

        return calcinfo

    def create_restart(self, force_restart=False, parent_folder_symlink=None):
        """
        Function to restart a calculation that was not completed before 
        (like max walltime reached...) i.e. not to restart a really FAILED calculation.
        Returns a calculation c2, with all links prepared but not stored in DB.
        To submit it, simply do::

          c2.store_all()
          c2.submit()

        .. deprecated:: 3.0
           Use the helper method :py:func:`aiida_quantumespresso.utils.restart.create_restart_neb` instead,
           that returns a calculation builder rather than a new, unstored calculation.


        :param bool force_restart: restart also if parent is not in FINISHED 
        state (e.g. FAILED, IMPORTED, etc.). Default=False.
        :param bool parent_folder_symlink: if True, symlinks are used
        instead of hard copies of the files. Default given by 
        self._default_symlink_usage.
        """
        from aiida_quantumespresso.utils.restart import clone_calculation
        import warnings
        warnings.warn('This method has been deprecated, use instead '
                      'aiida_quantumespresso.utils.restart.create_restart_neb()', DeprecationWarning)

        # Check the calculation's state using ``from_attribute=True`` to
        # correctly handle IMPORTED calculations.
        if not self.is_finished_ok:
            if not force_restart:
                raise InputValidationError(
                    "Calculation to be restarted must be finshed ok. Otherwise, use the force_restart flag")

        if parent_folder_symlink is None:
            parent_folder_symlink = self._default_symlink_usage

        inputs = self.get_incoming(link_type=LinkType.INPUT_CALC)

        old_inp_dict = inputs.get_node_by_label(self.get_linkname('neb_parameters')).get_dict()
        # add the restart flag
        old_inp_dict['PATH']['restart_mode'] = 'restart'
        inp_dict = Dict(dict=old_inp_dict)

        try:
            remote_folder = self.get_outgoing(node_class=RemoteData, link_label_filter='remote_folder').one().node
        except ValueError:
            raise InputValidationError("No or more than one output RemoteData found in calculation {}".format(self.pk))

        c2 = clone_calculation(self)

        #if not 'Restart' in c2.label:
        #    labelstring = c2.label + " Restart of {} {}.".format(
        #                                self.__class__.__name__,self.pk)
        #else:
        #    labelstring = " Restart of {} {}.".format(self.__class__.__name__,self.pk)
        #c2.label = labelstring.lstrip()

        # set the new links
        c2.use_neb_parameters(inp_dict)
        
        c2.use_pw_parameters(inputs.get_node_by_label(self.get_linkname('pw_parameters')))
        
        c2.use_first_structure(inputs.get_node_by_label(self.get_linkname('first_structure')))
        c2.use_last_structure(inputs.get_node_by_label(self.get_linkname('last_structure')))

        if self._use_kpoints:
            c2.use_kpoints(inputs.get_node_by_label(self.get_linkname('kpoints')))
        c2.use_code(inputs.get_node_by_label(self.get_linkname('code')))
        try:
            old_settings_dict = inputs.get_node_by_label(self.get_linkname('settings')).get_dict()
        except KeyError:
            old_settings_dict = {}
        if parent_folder_symlink is not None:
            old_settings_dict['PARENT_FOLDER_SYMLINK'] = parent_folder_symlink

        if old_settings_dict:  # if not empty dictionary
            settings = Dict(dict=old_settings_dict)
            c2.use_settings(settings)

        c2._set_parent_remotedata(remote_folder)
        # set links for pseudos
        for triple in self.get_incoming(node_class=UpfData).all():
            c2._add_link_from(triple.node, label=triple.link_label)

        # Add also the vdw table, if the parent had one
        try:
            old_vdw_table = inputs.get_node_by_label(self.get_linkname('vdw_table'))
        except KeyError:
            # No VdW table
            pass
        else:
            c2.use_vdw_table(old_vdw_table)

        return c2
