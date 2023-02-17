# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the pw.x code of Quantum ESPRESSO."""
import os
import warnings

from aiida import orm
from aiida.common.lang import classproperty
from aiida.plugins import factories

from aiida_quantumespresso.calculations import BasePwCpInputGenerator


class PwCalculation(BasePwCpInputGenerator):
    """`CalcJob` implementation for the pw.x code of Quantum ESPRESSO."""

    _automatic_namelists = {
        'scf': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
        'nscf': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
        'bands': ['CONTROL', 'SYSTEM', 'ELECTRONS'],
        'relax': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS'],
        'md': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS'],
        'vc-md': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'],
        'vc-relax': ['CONTROL', 'SYSTEM', 'ELECTRONS', 'IONS', 'CELL'],
    }

    # Keywords that cannot be set by the user but will be set by the plugin
    _blocked_keywords = [
        ('CONTROL', 'pseudo_dir'),
        ('CONTROL', 'outdir'),
        ('CONTROL', 'prefix'),
        ('SYSTEM', 'celldm'),
        ('SYSTEM', 'nat'),
        ('SYSTEM', 'ntyp'),
        ('SYSTEM', 'a'),
        ('SYSTEM', 'b'),
        ('SYSTEM', 'c'),
        ('SYSTEM', 'cosab'),
        ('SYSTEM', 'cosac'),
        ('SYSTEM', 'cosbc'),
    ]

    _use_kpoints = True

    # Not using symlink in pw to allow multiple nscf to run on top of the same scf
    _default_symlink_usage = False

    _ENABLED_PARALLELIZATION_FLAGS = ('npool', 'nband', 'ntg', 'ndiag')

    @classproperty
    def xml_filepaths(cls):
        """Return a list of XML output filepaths relative to the remote working directory that should be retrieved."""
        # pylint: disable=no-self-argument,not-an-iterable
        filepaths = []

        for filename in cls.xml_filenames:
            filepath = os.path.join(cls._OUTPUT_SUBFOLDER, f'{cls._PREFIX}.save', filename)
            filepaths.append(filepath)

        return filepaths

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('metadata.options.parser_name', valid_type=str, default='quantumespresso.pw')
        spec.input('metadata.options.without_xml', valid_type=bool, required=False, help='If set to `True` the parser '
            'will not fail if the XML file is missing in the retrieved folder.')
        spec.input('kpoints', valid_type=orm.KpointsData,
            help='kpoint mesh or kpoint path')
        spec.input('hubbard_file', valid_type=orm.SinglefileData, required=False,
            help='SinglefileData node containing the output Hubbard parameters from a HpCalculation')
        spec.inputs.validator = cls.validate_inputs

        spec.output('output_parameters', valid_type=orm.Dict,
            help='The `output_parameters` output node of the successful calculation.')
        spec.output('output_structure', valid_type=orm.StructureData, required=False,
            help='The `output_structure` output node of the successful calculation if present.')
        spec.output('output_trajectory', valid_type=orm.TrajectoryData, required=False)
        spec.output('output_band', valid_type=orm.BandsData, required=False,
            help='The `output_band` output node of the successful calculation if present.')
        spec.output('output_kpoints', valid_type=orm.KpointsData, required=False)
        spec.output('output_atomic_occupations', valid_type=orm.Dict, required=False)
        spec.default_output_node = 'output_parameters'

        # Unrecoverable errors: required retrieved files could not be read, parsed or are otherwise incomplete
        spec.exit_code(301, 'ERROR_NO_RETRIEVED_TEMPORARY_FOLDER',
            message='The retrieved temporary folder could not be accessed.')
        spec.exit_code(302, 'ERROR_OUTPUT_STDOUT_MISSING',
            message='The retrieved folder did not contain the required stdout output file.')
        spec.exit_code(303, 'ERROR_OUTPUT_XML_MISSING',
            message='The retrieved folder did not contain the required XML file.')
        spec.exit_code(304, 'ERROR_OUTPUT_XML_MULTIPLE',
            message='The retrieved folder contained multiple XML files.')
        spec.exit_code(305, 'ERROR_OUTPUT_FILES',
            message='Both the stdout and XML output files could not be read or parsed.')
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(311, 'ERROR_OUTPUT_STDOUT_PARSE',
            message='The stdout output file could not be parsed.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete probably because the calculation got interrupted.')
        spec.exit_code(320, 'ERROR_OUTPUT_XML_READ',
            message='The XML output file could not be read.')
        spec.exit_code(321, 'ERROR_OUTPUT_XML_PARSE',
            message='The XML output file could not be parsed.')
        spec.exit_code(322, 'ERROR_OUTPUT_XML_FORMAT',
            message='The XML output file has an unsupported format.')
        spec.exit_code(340, 'ERROR_OUT_OF_WALLTIME_INTERRUPTED',
            message='The calculation stopped prematurely because it ran out of walltime but the job was killed by the '
                    'scheduler before the files were safely written to disk for a potential restart.')
        spec.exit_code(350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION',
            message='The parser raised an unexpected exception: {exception}')

        # Significant errors but calculation can be used to restart
        spec.exit_code(400, 'ERROR_OUT_OF_WALLTIME',
            message='The calculation stopped prematurely because it ran out of walltime.')
        spec.exit_code(410, 'ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED',
            message='The electronic minimization cycle did not reach self-consistency.')

        spec.exit_code(461, 'ERROR_DEXX_IS_NEGATIVE',
            message='The code failed with negative dexx in the exchange calculation.')
        spec.exit_code(462, 'ERROR_COMPUTING_CHOLESKY',
            message='The code failed during the cholesky factorization.')
        spec.exit_code(463, 'ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED',
            message='Too many bands failed to converge during the diagonalization.')

        spec.exit_code(481, 'ERROR_NPOOLS_TOO_HIGH',
            message='The k-point parallelization "npools" is too high, some nodes have no k-points.')

        spec.exit_code(500, 'ERROR_IONIC_CONVERGENCE_NOT_REACHED',
            message='The ionic minimization cycle did not converge for the given thresholds.')
        spec.exit_code(501, 'ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF',
            message='Then ionic minimization cycle converged but the thresholds are exceeded in the final SCF.')
        spec.exit_code(502, 'ERROR_IONIC_CYCLE_EXCEEDED_NSTEP',
            message='The ionic minimization cycle did not converge after the maximum number of steps.')
        spec.exit_code(510, 'ERROR_IONIC_CYCLE_ELECTRONIC_CONVERGENCE_NOT_REACHED',
            message='The electronic minimization cycle failed during an ionic minimization cycle.')
        spec.exit_code(511, 'ERROR_IONIC_CONVERGENCE_REACHED_FINAL_SCF_FAILED',
            message='The ionic minimization cycle converged, but electronic convergence was not reached in the '
                    'final SCF.')
        spec.exit_code(520, 'ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE',
            message='The ionic minimization cycle terminated prematurely because of two consecutive failures in the '
                    'BFGS algorithm.')
        spec.exit_code(521, 'ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE',
            message='The ionic minimization cycle terminated prematurely because of two consecutive failures in the '
                    'BFGS algorithm and electronic convergence failed in the final SCF.')

        spec.exit_code(531, 'ERROR_CHARGE_IS_WRONG',
            message='The electronic minimization cycle did not reach self-consistency.')
        spec.exit_code(541, 'ERROR_SYMMETRY_NON_ORTHOGONAL_OPERATION',
            message='The variable cell optimization broke the symmetry of the k-points.')

        # Strong warnings about calculation results, but something tells us that you're ok with that
        spec.exit_code(710, 'WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED',
            message='The electronic minimization cycle did not reach self-consistency, but `scf_must_converge` '
                    'is `False` and/or `electron_maxstep` is 0.')
        # yapf: enable

    @staticmethod
    def validate_inputs_base(value, _):
        """Validate the top level namespace."""
        parameters = value['parameters'].get_dict()
        calculation_type = parameters.get('CONTROL', {}).get('calculation', 'scf')

        # Check that the restart input parameters are set correctly
        if calculation_type in ('nscf', 'bands'):
            if parameters.get('ELECTRONS', {}).get('startingpot', 'file') != 'file':
                return f'`startingpot` should be set to `file` for a `{calculation_type}` calculation.'
            if parameters.get('CONTROL', {}).get('restart_mode', 'from_scratch') != 'from_scratch':
                warnings.warn(f'`restart_mode` should be set to `from_scratch` for a `{calculation_type}` calculation.')
        elif 'parent_folder' in value:
            if not any([
                parameters.get('CONTROL', {}).get('restart_mode', None) == 'restart',
                parameters.get('ELECTRONS', {}).get('startingpot', None) == 'file',
                parameters.get('ELECTRONS', {}).get('startingwfc', None) == 'file'
            ]):
                warnings.warn(
                    '`parent_folder` input was provided for the `PwCalculation`, but no '
                    'input parameters are set to restart from these files.'
                )

    @classmethod
    def validate_inputs(cls, value, port_namespace):
        """Validate the top level namespace.

        Check that the restart input parameters are set correctly. In case of 'nscf' and 'bands' calculations, this
        means ``parent_folder`` is provided, ``startingpot`` is set to 'file' and ``restart_mode`` is 'from_scratch'.
        For other calculations, if the ``parent_folder`` is provided, the restart settings must be set to use some of
        the outputs.

        Note that the validator is split in two methods: ``validate_inputs`` and ``validate_inputs_base``. This is to
        facilitate work chains that wrap this calculation that will provide the ``parent_folder`` themselves and so do
        not require the user to provide it at launch of the work chain. This will fail because of the validation in this
        validator, however, which is why the rest of the logic is moved to ``validate_inputs_base``. The wrapping work
        chain can change the ``validate_input`` validator for ``validate_inputs_base`` thereby allowing the parent
        folder to be defined during the work chains runtime, while still keep the rest of the namespace validation.
        """
        result = super().validate_inputs(value, port_namespace)

        if result is not None:
            return result

        parameters = value['parameters'].get_dict()
        calculation_type = parameters.get('CONTROL', {}).get('calculation', 'scf')

        if calculation_type in ('nscf', 'bands'):
            if 'parent_folder' not in value:
                warnings.warn(
                    f'`parent_folder` not provided for `{calculation_type}` calculation. For work chains wrapping this '
                    'calculation, you can disable this warning by setting the validator of the `PwCalculation` port to '
                    '`PwCalculation.validate_inputs_base`.'
                )

        return cls.validate_inputs_base(value, port_namespace)

    @classproperty
    def filename_input_hubbard_parameters(cls):
        """Return the relative file name of the file containing the Hubbard parameters.

        .. note:: This only applies if they should be read from file instead of specified in the input file cards.
        .. warning:: Requires the aiida-quantumespresso-hp plugin to be installed
        """
        # pylint: disable=no-self-argument,no-self-use
        try:
            HpCalculation = factories.CalculationFactory('quantumespresso.hp')
        except Exception as exc:
            raise RuntimeError(
                'this is determined by the aiida-quantumespresso-hp plugin but it is not installed'
            ) from exc

        return HpCalculation.filename_input_hubbard_parameters

    @classmethod
    def input_helper(cls, *args, **kwargs):
        """Validate the provided keywords and prepare the inputs dictionary in a 'standardized' form.

        The standardization converts ints to floats when required, or if the flag `flat_mode` is specified,
        puts the keywords in the right namelists.

        This function calls :py:func:`aiida_quantumespresso.calculations.helpers.pw_input_helper`, see its docstring for
        further information.
        """
        from . import helpers
        return helpers.pw_input_helper(*args, **kwargs)
