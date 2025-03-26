# -*- coding: utf-8 -*-
"""Plugin to create a Quantum Espresso ph.x input file."""
import os
import warnings

from aiida import orm
from aiida.common import datastructures, exceptions
from aiida.common.warnings import AiidaDeprecationWarning
import numpy

from aiida_quantumespresso.calculations import _lowercase_dict, _uppercase_dict
from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry

from .base import CalcJob


class PhCalculation(CalcJob):
    """`CalcJob` implementation for the ph.x code of Quantum ESPRESSO."""

    # Keywords that cannot be set by the user but will be set by the plugin
    _blocked_keywords = [('INPUTPH', 'outdir'), ('INPUTPH', 'verbosity'), ('INPUTPH', 'prefix'), ('INPUTPH', 'fildyn'),
                         ('INPUTPH', 'ldisp'), ('INPUTPH', 'nq1'), ('INPUTPH', 'nq2'), ('INPUTPH', 'nq3'),
                         ('INPUTPH', 'qplot')]

    _use_kpoints = True

    _compulsory_namelists = ['INPUTPH']

    # Default input and output files
    _PREFIX = 'aiida'
    _DEFAULT_INPUT_FILE = 'aiida.in'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'
    _OUTPUT_XML_TENSOR_FILE_NAME = 'tensors.xml'
    _OUTPUT_SUBFOLDER = './out/'
    _FOLDER_DRHO = 'FILDRHO'
    _DRHO_PREFIX = 'drho'
    _DVSCF_PREFIX = 'dvscf'
    _DRHO_STAR_EXT = 'drho_rot'
    _FOLDER_DYNAMICAL_MATRIX = 'DYN_MAT'
    _VERBOSITY = 'high'
    _OUTPUT_DYNAMICAL_MATRIX_PREFIX = os.path.join(_FOLDER_DYNAMICAL_MATRIX, 'dynamical-matrix-')

    # Not using symlink in pw to allow multiple nscf to run on top of the same scf
    _default_symlink_usage = False

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('metadata.options.input_filename', valid_type=str, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=str, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.parser_name', valid_type=str, default='quantumespresso.ph')
        spec.input('metadata.options.withmpi', valid_type=bool, default=True)
        spec.input('qpoints', valid_type=orm.KpointsData, help='qpoint mesh')
        spec.input('parameters', valid_type=orm.Dict, help='')
        spec.input('settings', valid_type=orm.Dict, required=False, help='')
        spec.input('parent_folder', valid_type=orm.RemoteData,
            help='the folder of a completed `PwCalculation`')
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.default_output_node = 'output_parameters'

        # Unrecoverable errors: required retrieved files could not be read, parsed or are otherwise incomplete
        spec.exit_code(302, 'ERROR_OUTPUT_STDOUT_MISSING',
            message='The retrieved folder did not contain the required stdout output file.')
        spec.exit_code(305, 'ERROR_OUTPUT_FILES',
            message='Both the stdout and XML output files could not be read or parsed.')
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(311, 'ERROR_OUTPUT_STDOUT_PARSE',
            message='The stdout output file could not be parsed.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete probably because the calculation got interrupted.')
        spec.exit_code(350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION',
            message='The parser raised an unexpected exception: {exception}')
        spec.exit_code(360, 'ERROR_INCOMPATIBLE_FFT_GRID',
            message='The FFT grid is incompatible with the detected symmetries. Try using the lattice-specific '
                    '`ibrav` != 0 in the parent `pw.x` calculation.')
        spec.exit_code(361, 'ERROR_WRONG_REPRESENTATION',
            message=('The representation found seems to be wrong according to the detected symmetries. '
                     'Try using the lattice-specific `ibrav` != 0 in the parent `pw.x` calculation.'))

        # Significant errors but calculation can be used to restart
        spec.exit_code(400, 'ERROR_OUT_OF_WALLTIME',
            message='The calculation stopped prematurely because it ran out of walltime.')
        spec.exit_code(410, 'ERROR_CONVERGENCE_NOT_REACHED',
            message='The minimization cycle did not reach self-consistency.')
        spec.exit_code(462, 'ERROR_COMPUTING_CHOLESKY',
            message='The code failed during the cholesky factorization.')
        # yapf: enable

    def prepare_for_submission(self, folder):
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :py:class:`~aiida.common.datastructures.CalcInfo` instance.
        """
        # pylint: disable=too-many-statements,too-many-branches
        local_copy_list = []
        remote_copy_list = []
        remote_symlink_list = []

        if 'settings' in self.inputs:
            settings = _uppercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
            if 'ADDITIONAL_RETRIEVE_LIST' in settings:
                warnings.warn(
                    'The key `ADDITIONAL_RETRIEVE_LIST` in the settings input is deprecated and will be removed in '
                    'the future. Use the `CalcJob.metadata.options.additional_retrieve_list` input instead.',
                    AiidaDeprecationWarning
                )
        else:
            settings = {}

        parent_folder = self.inputs.parent_folder
        parent_calcs = parent_folder.base.links.get_incoming(node_class=orm.CalcJobNode).all()

        if not parent_calcs:
            raise exceptions.NotExistent(f'parent_folder<{parent_folder.pk}> has no parent calculation')
        elif len(parent_calcs) > 1:
            raise exceptions.UniquenessError(f'parent_folder<{parent_folder.pk}> has multiple parent calculations')

        parent_calc = parent_calcs[0].node

        # If the parent calculation is a `PhCalculation` we are restarting
        restart_flag = parent_calc.process_type == 'aiida.calculations:quantumespresso.ph'

        # Also, the parent calculation must be on the same computer
        if self.node.computer.uuid != parent_calc.computer.uuid:
            computer_label = parent_calc.computer.label
            raise exceptions.InputValidationError(
                f'Calculation has to be launched on the same computer as that of the parent: {computer_label}'
            )

        # put by default, default_parent_output_folder = ./out
        try:
            default_parent_output_folder = parent_calc.process_class._OUTPUT_SUBFOLDER  # pylint: disable=protected-access
        except AttributeError:
            try:
                default_parent_output_folder = parent_calc._get_output_folder()  # pylint: disable=protected-access
            except AttributeError as exception:
                msg = 'parent calculation does not have a default output subfolder'
                raise exceptions.InputValidationError(msg) from exception
        parent_calc_out_subfolder = settings.pop('PARENT_CALC_OUT_SUBFOLDER', default_parent_output_folder)

        # I put the first-level keys as uppercase (i.e., namelist and card names) and the second-level keys as lowercase
        parameters = _uppercase_dict(self.inputs.parameters.get_dict(), dict_name='parameters')
        parameters = {k: _lowercase_dict(v, dict_name=k) for k, v in parameters.items()}

        prepare_for_d3 = settings.pop('PREPARE_FOR_D3', False)
        if prepare_for_d3:
            self._blocked_keywords += [('INPUTPH', 'fildrho'), ('INPUTPH', 'drho_star%open'),
                                       ('INPUTPH', 'drho_star%ext'), ('INPUTPH', 'drho_star%dir')]

        for namelist, flag in self._blocked_keywords:
            if namelist in parameters:
                if flag in parameters[namelist]:
                    raise exceptions.InputValidationError(
                        f"Cannot specify explicitly the '{flag}' flag in the '{namelist}' namelist or card."
                    )

        if 'INPUTPH' not in parameters:
            raise exceptions.InputValidationError('required namelist INPUTPH not specified')

        parameters['INPUTPH']['outdir'] = self._OUTPUT_SUBFOLDER
        parameters['INPUTPH']['verbosity'] = self._VERBOSITY
        parameters['INPUTPH']['prefix'] = self._PREFIX
        parameters['INPUTPH']['fildyn'] = self._OUTPUT_DYNAMICAL_MATRIX_PREFIX

        prepare_for_epw = settings.pop('PREPARE_FOR_EPW', False)
        if prepare_for_epw:
            self._blocked_keywords += [('INPUTPH', 'fildvscf')]
            parameters['INPUTPH']['fildvscf'] = self._DVSCF_PREFIX

        if prepare_for_d3:
            parameters['INPUTPH']['fildrho'] = self._DRHO_PREFIX
            parameters['INPUTPH']['drho_star%open'] = True
            parameters['INPUTPH']['drho_star%ext'] = self._DRHO_STAR_EXT
            parameters['INPUTPH']['drho_star%dir'] = self._FOLDER_DRHO

        try:
            mesh, offset = self.inputs.qpoints.get_kpoints_mesh()

            if any(i != 0. for i in offset):
                raise NotImplementedError(
                    'Computation of phonons on a mesh with non zero offset is not implemented, at the level of ph.x'
                )

            parameters['INPUTPH']['ldisp'] = True
            parameters['INPUTPH']['nq1'] = mesh[0]
            parameters['INPUTPH']['nq2'] = mesh[1]
            parameters['INPUTPH']['nq3'] = mesh[2]

            postpend_text = None

        except AttributeError:
            # this is the case where no mesh was set. Maybe it's a list
            try:
                list_of_points = self.inputs.qpoints.get_kpoints(cartesian=True)
            except AttributeError as exception:
                # In this case, there are no info on the qpoints at all
                msg = 'Input `qpoints` contains neither a mesh nor a list of points'
                raise exceptions.InputValidationError(msg) from exception

            # change to 2pi/a coordinates
            lattice_parameter = numpy.linalg.norm(self.inputs.qpoints.cell[0])
            list_of_points *= lattice_parameter / (2. * numpy.pi)

            # add here the list of point coordinates
            if len(list_of_points) > 1:
                parameters['INPUTPH']['qplot'] = True
                parameters['INPUTPH']['ldisp'] = True
                postpend_text = f'{len(list_of_points)}\n'
                for points in list_of_points:
                    postpend_text += '{0:18.10f} {1:18.10f} {2:18.10f}  1\n'.format(*points)  # pylint: disable=consider-using-f-string

                # Note: the weight is fixed to 1, because ph.x calls these
                # things weights but they are not such. If they are going to
                # exist with the meaning of weights, they will be supported
            else:
                parameters['INPUTPH']['ldisp'] = False
                postpend_text = ''
                for points in list_of_points:
                    postpend_text += '{0:18.10f} {1:18.10f} {2:18.10f}\n'.format(*points)  # pylint: disable=consider-using-f-string

        # customized namelists, otherwise not present in the distributed ph code
        try:
            namelists_toprint = settings.pop('NAMELISTS')
            if not isinstance(namelists_toprint, list):
                raise exceptions.InputValidationError(
                    "The 'NAMELISTS' value, if specified in the settings input "
                    'node, must be a list of strings'
                )
        except KeyError:  # list of namelists not specified in the settings; do automatic detection
            namelists_toprint = self._compulsory_namelists

        # create a folder for the dynamical matrices
        if not restart_flag:  # if it is a restart, it will be copied over
            folder.get_subfolder(self._FOLDER_DYNAMICAL_MATRIX, create=True)

        with folder.open(self.metadata.options.input_filename, 'w') as infile:
            for namelist_name in namelists_toprint:
                infile.write(f'&{namelist_name}\n')
                # namelist content; set to {} if not present, so that we leave an empty namelist
                namelist = parameters.pop(namelist_name, {})
                for key, value in sorted(namelist.items()):
                    infile.write(convert_input_to_namelist_entry(key, value))
                infile.write('/\n')

            # add list of qpoints if required
            if postpend_text is not None:
                infile.write(postpend_text)

        if parameters:
            raise exceptions.InputValidationError(
                'The following namelists are specified in parameters, but are not valid namelists for the current type '
                f'of calculation: {",".join(list(parameters.keys()))}'
            )

        # copy the parent scratch
        symlink = settings.pop('PARENT_FOLDER_SYMLINK', self._default_symlink_usage)  # a boolean
        if symlink:
            # I create a symlink to each file/folder in the parent ./out
            folder.get_subfolder(self._OUTPUT_SUBFOLDER, create=True)

            remote_symlink_list.append((
                parent_folder.computer.uuid,
                os.path.join(parent_folder.get_remote_path(), parent_calc_out_subfolder, '*'), self._OUTPUT_SUBFOLDER
            ))

            # I also create a symlink for the ./pseudo folder
            # Remove this when the recover option of QE will be fixed (bug when trying to find pseudo file)
            remote_symlink_list.append((
                parent_folder.computer.uuid, os.path.join(parent_folder.get_remote_path(),
                                                          self._get_pseudo_folder()), self._get_pseudo_folder()
            ))
        else:
            # here I copy the whole folder ./out
            remote_copy_list.append((
                parent_folder.computer.uuid, os.path.join(parent_folder.get_remote_path(),
                                                          parent_calc_out_subfolder), self._OUTPUT_SUBFOLDER
            ))
            # I also copy the ./pseudo folder
            # Remove this when the recover option of QE will be fixed (bug when trying to find pseudo file)
            remote_copy_list.append((
                parent_folder.computer.uuid, os.path.join(parent_folder.get_remote_path(),
                                                          self._get_pseudo_folder()), self._get_pseudo_folder()
            ))

        if restart_flag:  # in this case, copy in addition also the dynamical matrices
            if symlink:
                remote_symlink_list.append((
                    parent_folder.computer.uuid,
                    os.path.join(parent_folder.get_remote_path(),
                                 self._FOLDER_DYNAMICAL_MATRIX), self._FOLDER_DYNAMICAL_MATRIX
                ))

            else:
                # copy the dynamical matrices
                # no need to copy the _ph0, since I copied already the whole ./out folder
                remote_copy_list.append((
                    parent_folder.computer.uuid,
                    os.path.join(parent_folder.get_remote_path(), self._FOLDER_DYNAMICAL_MATRIX), '.'
                ))

        # Create an `.EXIT` file if `only_initialization` flag in `settings` is set to `True`
        if settings.pop('ONLY_INITIALIZATION', False):
            with folder.open(f'{self._PREFIX}.EXIT', 'w') as handle:
                handle.write('\n')

        codeinfo = datastructures.CodeInfo()
        codeinfo.cmdline_params = (list(settings.pop('CMDLINE', [])) + ['-in', self.metadata.options.input_filename])
        codeinfo.stdout_name = self.metadata.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid

        calcinfo = datastructures.CalcInfo()
        calcinfo.uuid = str(self.uuid)
        calcinfo.codes_info = [codeinfo]
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.remote_symlink_list = remote_symlink_list

        # Retrieve by default the output file and the xml file
        filepath_xml_tensor = os.path.join(self._OUTPUT_SUBFOLDER, '_ph0', f'{self._PREFIX}.phsave')
        calcinfo.retrieve_list = []
        calcinfo.retrieve_list.append(self.metadata.options.output_filename)
        calcinfo.retrieve_list.append(self._FOLDER_DYNAMICAL_MATRIX)
        calcinfo.retrieve_list.append(os.path.join(filepath_xml_tensor, self._OUTPUT_XML_TENSOR_FILE_NAME))
        calcinfo.retrieve_list += settings.pop('ADDITIONAL_RETRIEVE_LIST', [])

        if settings:
            unknown_keys = ', '.join(list(settings.keys()))
            raise exceptions.InputValidationError(f'`settings` contained unexpected keys: {unknown_keys}')

        return calcinfo

    @staticmethod
    def _get_pseudo_folder():
        """Get the calculation-specific pseudo folder (relative path).

        Default given by PwCalculation._PSEUDO_SUBFOLDER
        """
        return PwCalculation._PSEUDO_SUBFOLDER  # pylint: disable=protected-access
