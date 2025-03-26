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
        # pylint: disable=too-many-statements
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
        spec.exit_code(303, 'ERROR_OUTPUT_XML_MISSING',
            message='The retrieved folder did not contain the required XML file.')
        spec.exit_code(304, 'ERROR_OUTPUT_XML_MULTIPLE',
            message='The retrieved folder contained multiple XML files.')
        spec.exit_code(305, 'ERROR_OUTPUT_FILES',
            message='Both the stdout and XML output files could not be read or parsed.')
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
        spec.exit_code(360, 'ERROR_G_PAR',
            message='The code failed in finding a valid reciprocal lattice vector.')

        # Significant errors but calculation can be used to restart
        spec.exit_code(410, 'ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED',
            message='The electronic minimization cycle did not reach self-consistency.')

        spec.exit_code(461, 'ERROR_DEXX_IS_NEGATIVE',
            message='The code failed with negative dexx in the exchange calculation.')
        spec.exit_code(462, 'ERROR_COMPUTING_CHOLESKY',
            message='The code failed during the cholesky factorization.')
        spec.exit_code(463, 'ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED',
            message='Too many bands failed to converge during the diagonalization.')
        spec.exit_code(464, 'ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE',
            message='The S matrix was found to be not positive definite.')
        spec.exit_code(465, 'ERROR_ZHEGVD_FAILED',
            message='The `zhegvd` failed in the PPCG diagonalization.')
        spec.exit_code(466, 'ERROR_QR_FAILED',
            message='The `[Q, R] = qr(X, 0)` failed in the PPCG diagonalization.')
        spec.exit_code(467, 'ERROR_EIGENVECTOR_CONVERGENCE',
            message='The eigenvector failed to converge.')
        spec.exit_code(468, 'ERROR_BROYDEN_FACTORIZATION',
            message='The factorization in the Broyden routine failed.')

        spec.exit_code(481, 'ERROR_NPOOLS_TOO_HIGH',
            message='The k-point parallelization "npools" is too high, some nodes have no k-points.')

        spec.exit_code(500, 'ERROR_IONIC_CONVERGENCE_NOT_REACHED',
            message='The ionic minimization cycle did not converge for the given thresholds.')
        spec.exit_code(501, 'ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF',
            message='Then ionic minimization cycle converged but the thresholds are exceeded in the final SCF.')
        spec.exit_code(502, 'ERROR_IONIC_CYCLE_EXCEEDED_NSTEP',
            message='The ionic minimization cycle did not converge after the maximum number of steps.')
        spec.exit_code(503, 'ERROR_IONIC_INTERRUPTED_PARTIAL_TRAJECTORY',
            message='The ionic minimization cycle did not finish because the calculation was interrupted but a partial '
                    'trajectory and output structure was successfully parsed which can be used for a restart.')
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
            message='The difference between the total charge and the number of electrons exceeds the threshold. '
                    'Smearing might be required, check the output file for details.'
            )
        spec.exit_code(541, 'ERROR_SYMMETRY_NON_ORTHOGONAL_OPERATION',
            message='The variable cell optimization broke the symmetry of the k-points.')
        spec.exit_code(542, 'ERROR_RADIAL_FFT_SIGNIFICANT_VOLUME_CONTRACTION',
            message=('The cell relaxation caused a significant volume contraction '
                     'and there is not enough space allocated for radial FFT.'))

        # Strong warnings about calculation results, but something tells us that you're ok with that
        spec.exit_code(710, 'WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED',
            message='The electronic minimization cycle did not reach self-consistency, but `scf_must_converge` '
                    'is `False` and/or `electron_maxstep` is 0.')
        # yapf: enable

    @classmethod
    def validate_inputs(cls, value, port_namespace):
        """Validate the top level namespace.

        Check that the restart input parameters are set correctly. In case of 'nscf' and 'bands' calculations, this
        means ``parent_folder`` is provided. For other calculations, if the ``parent_folder`` is provided, the restart
        settings must be set to use some of the outputs.

        Note that the validator will only check the logic in case the ``parent_folder`` is a port in the
        ``port_namespace``. This is because the ``PwCalculation`` can be wrapped inside a work chain that only provides
        the ``parent_folder`` input at a later step in the outline. To avoid raising any warnings, such a work chain
        must exclude the ``parent_folder`` port when exposing the inputs of the ``PwCalculation``.
        """
        from aiida.engine.processes.calcjobs.calcjob import validate_calc_job

        result = super().validate_inputs(value, port_namespace)

        if result is not None:
            return result

        result = validate_calc_job(value, port_namespace)

        if result is not None:
            return result

        parameters = value['parameters'].get_dict()
        calculation_type = parameters.get('CONTROL', {}).get('calculation', 'scf')

        if 'parent_folder' not in port_namespace:
            return

        # If a `parent_folder` input is provided, make sure the inputs are set to restart
        if 'parent_folder' in value and calculation_type not in ('nscf', 'bands'):
            if not any([
                parameters.get('CONTROL', {}).get('restart_mode', None) == 'restart',
                parameters.get('ELECTRONS', {}).get('startingpot', None) == 'file',
                parameters.get('ELECTRONS', {}).get('startingwfc', None) == 'file'
            ]):
                warnings.warn(
                    f'`parent_folder` input was provided for the `{calculation_type}` `PwCalculation`, but no input '
                    'parameters were provided to restart from this folder.\n\n'
                    'Please set one of the following in the input parameters:\n'
                    "    parameters['CONTROL']['restart_mode'] = 'restart'\n"
                    "    parameters['ELECTRONS']['startingpot'] = 'file'\n"
                    "    parameters['ELECTRONS']['startingwfc'] = 'file'\n", UserWarning
                )

        if calculation_type in ('nscf', 'bands'):
            if 'parent_folder' not in value:
                warnings.warn(
                    f'`parent_folder` not provided for `{calculation_type}` calculation. For work chains wrapping this '
                    'calculation, you can disable this warning by excluding the `parent_folder` when exposing the '
                    'inputs of the `PwCalculation`.', UserWarning
                )

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
