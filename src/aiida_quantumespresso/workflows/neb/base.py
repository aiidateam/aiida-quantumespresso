# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO neb.x calculation with automated error handling and restarts."""
from aiida import orm
from aiida.common import AttributeDict, InputValidationError
from aiida.engine import BaseRestartWorkChain, ProcessHandlerReport, process_handler, while_

from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import create_kpoints_from_distance
from aiida_quantumespresso.common.types import ElectronicType, RestartType, SpinType
from aiida_quantumespresso.utils.defaults.calculation import pw as qe_defaults

from ...calculations.neb import NebCalculation
from ...workflows.pw.base import PwBaseWorkChain
from ..protocols.utils import ProtocolMixin


class NebBaseWorkChain(ProtocolMixin, BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO neb.x calculation with automated error handling and restarts."""

    # pylint: disable=too-many-public-methods, too-many-statements

    _process_class = NebCalculation

    defaults = AttributeDict({
        'qe': qe_defaults,
        'delta_factor_mixing_beta': 0.8,
    })

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.expose_inputs(NebCalculation, namespace='neb', exclude=('pw.kpoints', 'first_structure', 'last_structure'))
        spec.input('kpoints', valid_type=orm.KpointsData, required=False,
            help='An explicit k-points list or mesh. Either this or `kpoints_distance` has to be provided.')
        spec.input('kpoints_distance', valid_type=orm.Float, required=False,
            help='The minimum desired distance in 1/Å between k-points in reciprocal space. The explicit k-points will '
                 'be generated automatically by a calculation function based on the input structure.')
        spec.input('kpoints_force_parity', valid_type=orm.Bool, required=False,
            help='Optional input when constructing the k-points based on a desired `kpoints_distance`. Setting this to '
                 '`True` will force the k-point mesh to have an even number of points along each lattice vector except '
                 'for any non-periodic directions.')

        spec.outline(
            cls.setup,
            cls.validate_kpoints,
            while_(cls.should_run_process)(
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )

        spec.expose_outputs(NebCalculation)

        spec.exit_code(202, 'ERROR_INVALID_INPUT_KPOINTS',
            message='Neither the `kpoints` nor the `kpoints_distance` input was specified.')

        spec.exit_code(300, 'ERROR_UNRECOVERABLE_FAILURE',
            message='The calculation failed with an unidentified unrecoverable error.')
        spec.exit_code(310, 'ERROR_KNOWN_UNRECOVERABLE_FAILURE',
            message='The calculation failed with a known unrecoverable error.')
        # yapf: enable

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        from ..protocols import pw as pw_protocols
        return files(pw_protocols) / 'base.yaml'

    @classmethod
    def get_builder_from_protocol(
        cls,
        code,
        images,
        protocol=None,
        overrides=None,
        electronic_type=ElectronicType.METAL,
        spin_type=SpinType.NONE,
        initial_magnetic_moments=None,
        options=None,
        **kargs
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param images: the ``TrajectoryData`` instance to use for initial guess of NEB images.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param electronic_type: indicate the electronic character of the system through ``ElectronicType`` instance.
        :param spin_type: indicate the spin polarization type to use through a ``SpinType`` instance.
        :param initial_magnetic_moments: optional dictionary that maps the initial magnetic moment of each kind to a
            desired value for a spin polarized calculation. Note that in case the ``starting_magnetization`` is also
            provided in the ``overrides``, this takes precedence over the values provided here. In case neither is
            provided and ``spin_type == SpinType.COLLINEAR``, an initial guess for the magnetic moments is used.
        :param options: A dictionary of options that will be recursively set for the ``metadata.options`` input of all
            the ``CalcJobs`` that are nested in this work chain.
        :return: a process builder instance with all inputs defined ready for launch.
        """

        pw_base = PwBaseWorkChain.get_builder_from_protocol(
            code,
            images.get_step_structure(-1),
            protocol=protocol,
            overrides=overrides,
            electronic_type=electronic_type,
            spin_type=spin_type,
            initial_magnetic_moments=initial_magnetic_moments,
            options=options,
            **kargs
        )
        #pylint: disable=no-member
        builder = cls.get_builder()
        builder.neb.code = code
        builder.neb.images = images
        builder.neb.pw.pseudos = pw_base.pw.pseudos
        builder.neb.pw.parameters = pw_base.pw.parameters
        builder.neb.metadata.options = pw_base.pw.metadata.options

        if 'kpoints' in pw_base:
            builder.kpoints = pw_base['kpoints']
        else:
            builder.kpoints_distance = orm.Float(pw_base['kpoints_distance'])
        builder.kpoints_force_parity = orm.Bool(pw_base['kpoints_force_parity'])
        builder.max_iterations = orm.Int(pw_base['max_iterations'])
        # pylint: enable=no-member
        return builder

    def setup(self):
        """Call the ``setup`` of the ``BaseRestartWorkChain`` and create the inputs dictionary in ``self.ctx.inputs``.

        This ``self.ctx.inputs`` dictionary will be used by the ``BaseRestartWorkChain`` to submit the calculations
        in the internal loop.

        The ``parameters`` and ``settings`` input ``Dict`` nodes are converted into a regular dictionary and the
        default namelists for the ``parameters`` are set to empty dictionaries if not specified.
        """
        super().setup()

        if 'images' not in self.inputs.neb:
            raise InputValidationError('Input `images` are required.')

        self.ctx.inputs = AttributeDict(self.exposed_inputs(NebCalculation, 'neb'))

        self.ctx.inputs.parameters = self.ctx.inputs.parameters.get_dict()
        self.ctx.inputs.pw.parameters = self.ctx.inputs.pw.parameters.get_dict()
        self.ctx.inputs.pw.parameters.setdefault('CONTROL', {})
        self.ctx.inputs.pw.parameters.setdefault('ELECTRONS', {})
        self.ctx.inputs.pw.parameters.setdefault('SYSTEM', {})

        self.ctx.inputs.pw.parameters['CONTROL'].setdefault('calculation', 'scf')

    def validate_kpoints(self):
        """Validate the inputs related to k-points.

        Either an explicit `KpointsData` with given mesh/path, or a desired k-points distance should be specified. In
        the case of the latter, the `KpointsData` will be constructed for the input `StructureData` using the
        `create_kpoints_from_distance` calculation function.
        """
        if all(key not in self.inputs for key in ['kpoints', 'kpoints_distance']):
            return self.exit_codes.ERROR_INVALID_INPUT_KPOINTS

        try:
            kpoints = self.inputs.kpoints
        except AttributeError:
            inputs = {
                'structure': self.inputs.neb.images.get_step_structure(-1),
                'distance': self.inputs.kpoints_distance,
                'force_parity': self.inputs.get('kpoints_force_parity', orm.Bool(False)),
                'metadata': {
                    'call_link_label': 'create_kpoints_from_distance'
                }
            }
            kpoints = create_kpoints_from_distance(**inputs)  # pylint: disable=unexpected-keyword-arg

        self.ctx.inputs.pw.kpoints = kpoints

    def set_restart_type(self, restart_type, parent_folder=None):
        """Set the restart type for the next iteration."""

        if parent_folder is None and restart_type != RestartType.FROM_SCRATCH:
            raise ValueError('When not restarting from scratch, a `parent_folder` must be provided.')

        if restart_type == RestartType.FROM_SCRATCH:
            self.ctx.inputs.parameters['PATH']['restart_mode'] = 'from_scratch'
            self.ctx.inputs.pop('parent_folder', None)

        elif restart_type == RestartType.FULL:
            self.ctx.inputs.parameters['PATH']['restart_mode'] = 'restart'
            self.ctx.inputs.parent_folder = parent_folder

    def report_error_handled(self, calculation, action):
        """Report an action taken for a calculation that has failed.

        This should be called in a registered error handler if its condition is met and an action was taken.

        :param calculation: the failed calculation node
        :param action: a string message with the action taken
        """
        arguments = [calculation.process_label, calculation.pk, calculation.exit_status, calculation.exit_message]
        self.report('{}<{}> failed with exit status {}: {}'.format(*arguments))
        self.report(f'Action taken: {action}')

    @process_handler(
        priority=585,
        exit_codes=[
            NebCalculation.exit_codes.ERROR_COMPUTING_CHOLESKY,
            NebCalculation.exit_codes.ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED,
            NebCalculation.exit_codes.ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE,
            NebCalculation.exit_codes.ERROR_ZHEGVD_FAILED,
            NebCalculation.exit_codes.ERROR_QR_FAILED,
            NebCalculation.exit_codes.ERROR_EIGENVECTOR_CONVERGENCE,
            NebCalculation.exit_codes.ERROR_BROYDEN_FACTORIZATION,
        ]
    )
    def handle_diagonalization_errors(self, calculation):
        """Handle known issues related to the diagonalization.

        We use the following strategy. When a diagonalization algorithm fails, we try using an other one
        still not used. Conjugate gradient (CG) is kept as last option, as it is the slowest among the
        available ones, but on the contrary it is the most stable as well, thus kept as `last resort`.

        Once the error handler has tried all ``diagonalization`` options, abort.
        """
        current = self.ctx.inputs.pw.parameters['ELECTRONS'].get('diagonalization', 'david')

        if 'diagonalizations' not in self.ctx:
            # Initialize a list to track diagonalisations that haven't been tried in reverse order or preference
            self.ctx.diagonalizations = [value for value in ['cg', 'paro', 'ppcg', 'david'] if value != current.lower()]

        try:
            new = self.ctx.diagonalizations.pop()
            self.ctx.inputs.pw.parameters['ELECTRONS']['diagonalization'] = new
            action = f'found diagonalization issues for ``{current}``, switching to ``{new}`` diagonalization.'
            self.report_error_handled(calculation, action)
            return ProcessHandlerReport(True)
        except IndexError:
            action = 'found diagonalization issues but already exploited all supported algorithms, aborting...'
            self.report_error_handled(calculation, action)
            return ProcessHandlerReport(True, self.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE)

    @process_handler(priority=575, exit_codes=[
        NebCalculation.exit_codes.ERROR_NEB_INTERRUPTED_PARTIAL_TRAJECTORY,
    ])
    def handle_neb_interrupted_partial_trajectory(self, calculation):
        """Handle `ERROR_NEB_INTERRUPTED_PARTIAL_TRAJECTORY` and exit code.

        In this case the calculation shut down cleanly and we can do a full restart.
        """

        self.set_restart_type(RestartType.FULL, calculation.outputs.remote_folder)
        self.report_error_handled(calculation, "restarting in full with `CONTROL.restart_mode` = 'restart'")

        return ProcessHandlerReport(True)

    @process_handler(priority=575, exit_codes=[
        NebCalculation.exit_codes.ERROR_NEB_CYCLE_EXCEEDED_NSTEP,
    ])
    def handle_neb_cycle_exceeded_nstep(self, calculation):
        """Handle `ERROR_NEB_CYCLE_EXCEEDED_NSTEP` and exit code.

        In this case the calculation shut down cleanly and we can do a full restart.
        """
        self.ctx.inputs.parameters['PATH'].setdefault('nstep_path', 1)
        input_nsteps = self.inputs.neb.parameters['PATH']['nstep_path'] if 'nstep_path' in self.inputs.neb.parameters[
            'PATH'] else 1
        self.ctx.inputs.parameters['PATH']['nstep_path'] += input_nsteps

        self.set_restart_type(RestartType.FULL, calculation.outputs.remote_folder)
        self.report_error_handled(calculation, "restarting in full with `CONTROL.restart_mode` = 'restart'")

        return ProcessHandlerReport(True)

    @process_handler(priority=410, exit_codes=[
        NebCalculation.exit_codes.ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED,
    ])
    def handle_electronic_convergence_not_reached(self, calculation):
        """Handle `ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED` error.

        Decrease the mixing beta and fully restart from the previous calculation.
        """
        factor = self.defaults.delta_factor_mixing_beta
        mixing_beta = self.ctx.inputs.pw.parameters.get('ELECTRONS',
                                                        {}).get('mixing_beta', self.defaults.qe.mixing_beta)
        mixing_beta_new = mixing_beta * factor

        self.ctx.inputs.pw.parameters['ELECTRONS']['mixing_beta'] = mixing_beta_new
        action = f'reduced beta mixing from {mixing_beta} to {mixing_beta_new} and restarting from the last calculation'

        self.set_restart_type(RestartType.FULL, calculation.outputs.remote_folder)
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True)

    @process_handler(
        priority=400, exit_codes=[
            NebCalculation.exit_codes.ERROR_NEB_INTERRUPTED_WITHOUT_PARTIAL_TRAJECTORY,
        ]
    )
    def handle_neb_interrupted_without_partial_trajectory(self, calculation):
        """Handle `ERROR_NEB_INTERRUPTED_WITHOUT_PARTIAL_TRAJECTORY` error.

        In this case the calculation was interrupted before completing the first NEB minimization step,
        so we cannot retrieve any partial trajectory, no way to restart the calculation.
        Probably the walltime was too short.
        """
        action = 'Calculation was interrupted before completing the first NEB minimization step.'
        action += 'An increase of the walltime is probably needed. Aborting...'
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True, self.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE)
