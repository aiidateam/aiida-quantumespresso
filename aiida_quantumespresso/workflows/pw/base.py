# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO pw.x calculation with automated error handling and restarts."""
from __future__ import absolute_import

from aiida import orm
from aiida.common import AttributeDict, exceptions
from aiida.engine import ToContext, if_, while_
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.common.workchain.utils import register_error_handler, ErrorHandlerReport
from aiida_quantumespresso.common.workchain.base.restart import BaseRestartWorkChain
from aiida_quantumespresso.utils.defaults.calculation import pw as qe_defaults
from aiida_quantumespresso.utils.mapping import update_mapping, prepare_process_inputs
from aiida_quantumespresso.utils.pseudopotential import validate_and_prepare_pseudos_inputs
from aiida_quantumespresso.utils.resources import get_default_options, get_pw_parallelization_parameters
from aiida_quantumespresso.utils.resources import cmdline_remove_npools, create_scheduler_resources
from aiida_quantumespresso.workflows.functions.create_kpoints_from_distance import create_kpoints_from_distance

PwCalculation = CalculationFactory('quantumespresso.pw')


class PwBaseWorkChain(BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO pw.x calculation with automated error handling and restarts."""

    _calculation_class = PwCalculation
    _error_handler_entry_point = 'aiida_quantumespresso.workflow_error_handlers.pw.base'

    defaults = AttributeDict({
        'qe': qe_defaults,
        'delta_threshold_degauss': 30,
        'delta_factor_degauss': 0.1,
        'delta_factor_mixing_beta': 0.8,
        'delta_factor_max_seconds': 0.95,
        'delta_factor_nbnd': 0.05,
        'delta_minimum_nbnd': 4,
    })

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super(PwBaseWorkChain, cls).define(spec)
        spec.expose_inputs(PwCalculation, namespace='pw', exclude=('kpoints',))
        spec.input('pw.metadata.options.resources', valid_type=dict, required=False)
        spec.input('kpoints', valid_type=orm.KpointsData, required=False,
            help='An explicit k-points list or mesh. Either this or `kpoints_distance` has to be provided.')
        spec.input('kpoints_distance', valid_type=orm.Float, required=False,
            help='The minimum desired distance in 1/Å between k-points in reciprocal space. The explicit k-points will '
                 'be generated automatically by a calculation function based on the input structure.')
        spec.input('kpoints_force_parity', valid_type=orm.Bool, required=False,
            help='Optional input when constructing the k-points based on a desired `kpoints_distance`. Setting this to '
                 '`True` will force the k-point mesh to have an even number of points along each lattice vector except '
                 'for any non-periodic directions.')
        spec.input('pseudo_family', valid_type=orm.Str, required=False,
            help='An alternative to specifying the pseudo potentials manually in `pseudos`: one can specify the name '
                 'of an existing pseudo potential family and the work chain will generate the pseudos automatically '
                 'based on the input structure.')
        spec.input('automatic_parallelization', valid_type=orm.Dict, required=False,
            help='When defined, the work chain will first launch an initialization calculation to determine the '
                 'dimensions of the problem, and based on this it will try to set optimal parallelization flags.')

        spec.outline(
            cls.setup,
            cls.validate_parameters,
            cls.validate_kpoints,
            cls.validate_pseudos,
            cls.validate_resources,
            if_(cls.should_run_init)(
                cls.validate_init_inputs,
                cls.run_init,
                cls.inspect_init,
            ),
            while_(cls.should_run_calculation)(
                cls.prepare_calculation,
                cls.run_calculation,
                cls.inspect_calculation,
            ),
            cls.results,
        )

        spec.expose_outputs(PwCalculation)
        spec.output('automatic_parallelization', valid_type=orm.Dict, required=False,
            help='The results of the automatic parallelization analysis if performed.')

        spec.exit_code(201, 'ERROR_INVALID_INPUT_PSEUDO_POTENTIALS',
            message='The explicit `pseudos` or `pseudo_family` could not be used to get the necessary pseudos.')
        spec.exit_code(202, 'ERROR_INVALID_INPUT_KPOINTS',
            message='Neither the `kpoints` nor the `kpoints_distance` input was specified.')
        spec.exit_code(203, 'ERROR_INVALID_INPUT_RESOURCES',
            message='Neither the `options` nor `automatic_parallelization` input was specified.')
        spec.exit_code(204, 'ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED',
            message='The `metadata.options` did not specify both `resources.num_machines` and `max_wallclock_seconds`.')
        spec.exit_code(210, 'ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_MISSING_KEY',
            message='Required key for `automatic_parallelization` was not specified.')
        spec.exit_code(211, 'ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_UNRECOGNIZED_KEY',
            message='Unrecognized keys were specified for `automatic_parallelization`.')
        spec.exit_code(300, 'ERROR_UNRECOVERABLE_FAILURE',
            message='The calculation failed with an unrecoverable error.')
        spec.exit_code(320, 'ERROR_INITIALIZATION_CALCULATION_FAILED',
            message='The initialization calculation failed.')
        spec.exit_code(501, 'ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF',
            message='Then ionic minimization cycle converged but the thresholds are exceeded in the final SCF.')

    def setup(self):
        """Call the `setup` of the `BaseRestartWorkChain` and then create the inputs dictionary in `self.ctx.inputs`.

        This `self.ctx.inputs` dictionary will be used by the `BaseRestartWorkChain` to submit the calculations in the
        internal loop.
        """
        super(PwBaseWorkChain, self).setup()
        self.ctx.inputs = AttributeDict(self.exposed_inputs(PwCalculation, 'pw'))

    def validate_parameters(self):
        """Validate inputs that might depend on each other and cannot be validated by the spec.

        Also define dictionary `inputs` in the context, that will contain the inputs for the calculation that will be
        launched in the `run_calculation` step.
        """
        self.ctx.inputs.parameters = self.ctx.inputs.parameters.get_dict()
        self.ctx.inputs.settings = self.ctx.inputs.settings.get_dict() if 'settings' in self.ctx.inputs else {}

        if 'parent_folder' in self.ctx.inputs:
            self.ctx.restart_calc = self.ctx.inputs.parent_folder.creator

        self.ctx.inputs.parameters.setdefault('CONTROL', {})
        self.ctx.inputs.parameters['CONTROL'].setdefault('calculation', 'scf')

    def validate_kpoints(self):
        """Validate the inputs related to k-points.

        Either an explicit `KpointsData` with given mesh/path, or a desired k-points distance should be specified. In
        the case of the latter, the `KpointsData` will be constructed for the input `StructureData` using the
        `create_kpoints_from_distance` calculation function.
        """
        if all([key not in self.inputs for key in ['kpoints', 'kpoints_distance']]):
            return self.exit_codes.ERROR_INVALID_INPUT_KPOINTS

        try:
            kpoints = self.inputs.kpoints
        except AttributeError:
            inputs = {
                'structure': self.inputs.pw.structure,
                'distance': self.inputs.kpoints_distance,
                'force_parity': self.inputs.get('kpoints_force_parity', orm.Bool(False)),
                'metadata': {'call_link_label': 'create_kpoints_from_distance'}
            }
            kpoints = create_kpoints_from_distance(**inputs)  # pylint: disable=unexpected-keyword-arg

        self.ctx.inputs.kpoints = kpoints

    def validate_pseudos(self):
        """Validate the inputs related to pseudopotentials.

        Either the pseudo potentials should be defined explicitly in the `pseudos` namespace, or alternatively, a family
        can be defined in `pseudo_family` that will be used together with the input `StructureData` to generate the
        required mapping.
        """
        structure = self.ctx.inputs.structure
        pseudos = self.ctx.inputs.get('pseudos', None)
        pseudo_family = self.inputs.get('pseudo_family', None)

        try:
            self.ctx.inputs.pseudos = validate_and_prepare_pseudos_inputs(structure, pseudos, pseudo_family)
        except ValueError as exception:
            self.report('{}'.format(exception))
            return self.exit_codes.ERROR_INVALID_INPUT_PSEUDO_POTENTIALS

    def validate_resources(self):
        """Validate the inputs related to the resources.

        One can omit the normally required `options.resources` input for the `PwCalculation`, as long as the input
        `automatic_parallelization` is specified. If this is not the case, the `metadata.options` should at least
        contain the options `resources` and `max_wallclock_seconds`, where `resources` should define the `num_machines`.
        """
        if 'automatic_parallelization' not in self.inputs and 'options' not in self.ctx.inputs.metadata:
            return self.exit_codes.ERROR_INVALID_INPUT_RESOURCES

        # If automatic parallelization is not enabled, we better make sure that the options satisfy minimum requirements
        if 'automatic_parallelization' not in self.inputs:
            num_machines = self.ctx.inputs.metadata.options.get('resources', {}).get('num_machines', None)
            max_wallclock_seconds = self.ctx.inputs.metadata.options.get('max_wallclock_seconds', None)

            if num_machines is None or max_wallclock_seconds is None:
                return self.exit_codes.ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED

            self.set_max_seconds(max_wallclock_seconds)

    def set_max_seconds(self, max_wallclock_seconds):
        """Set the `max_seconds` to a fraction of `max_wallclock_seconds` option to prevent out-of-walltime problems.

        :param max_wallclock_seconds: the maximum wallclock time that will be set in the scheduler settings.
        """
        max_seconds_factor = self.defaults.delta_factor_max_seconds
        max_seconds = max_wallclock_seconds * max_seconds_factor
        self.ctx.inputs.parameters['CONTROL']['max_seconds'] = max_seconds

    def should_run_init(self):
        """Return whether an initialization calculation should be run.

        :return: boolean, `True` if `automatic_parallelization` was specified in the inputs, `False` otherwise.
        """
        return 'automatic_parallelization' in self.inputs

    def validate_init_inputs(self):
        """Validate the inputs that are required for the initialization calculation.

        The `automatic_parallelization` input expects a `Dict` node with the following keys:

            * max_wallclock_seconds
            * target_time_seconds
            * max_num_machines

        If any of these keys are not set or any superfluous keys are specified, the workchain will abort.
        """
        parallelization = self.inputs.automatic_parallelization.get_dict()

        expected_keys = ['max_wallclock_seconds', 'target_time_seconds', 'max_num_machines']
        received_keys = [(key, parallelization.get(key, None)) for key in expected_keys]
        remaining_keys = [key for key in parallelization.keys() if key not in expected_keys]

        for key, value in [(key, value) for key, value in received_keys if value is None]:
            self.report('required key "{}" in automatic_parallelization input not found'.format(key))
            return self.exit_codes.ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_MISSING_KEY

        if remaining_keys:
            self.report('detected unrecognized keys in the automatic_parallelization input: {}'.format(
                ' '.join(remaining_keys)))
            return self.exit_codes.ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_UNRECOGNIZED_KEY

        # Add the calculation mode to the automatic parallelization dictionary
        self.ctx.automatic_parallelization = {
            'max_wallclock_seconds': parallelization['max_wallclock_seconds'],
            'target_time_seconds': parallelization['target_time_seconds'],
            'max_num_machines': parallelization['max_num_machines'],
            'calculation_mode': self.ctx.inputs.parameters['CONTROL']['calculation']
        }

        options = self.ctx.inputs.metadata['options']
        options.setdefault('resources', {})['num_machines'] = parallelization['max_num_machines']
        options['max_wallclock_seconds'] = parallelization['max_wallclock_seconds']

    def run_init(self):
        """Run an initialization `PwCalculation` that will exit after the preamble.

        In the preamble, all the relevant dimensions of the problem are computed which allows us to make an estimate of
        the required resources and what parallelization flags need to be set.
        """
        inputs = self.ctx.inputs

        # Set the initialization flag and the initial default options
        inputs.settings['ONLY_INITIALIZATION'] = True
        inputs.metadata['options'] = update_mapping(inputs.metadata['options'], get_default_options())

        # Prepare the final input dictionary
        inputs = prepare_process_inputs(PwCalculation, inputs)
        running = self.submit(PwCalculation, **inputs)

        self.report('launching initialization PwCalculation<{}>'.format(running.pk))

        return ToContext(calculation_init=running)

    def inspect_init(self):
        """Use the initialization `PwCalculation` to determine the required resource and parallelization settings."""
        calculation = self.ctx.calculation_init

        if not calculation.is_finished_ok:
            return self.exit_codes.ERROR_INITIALIZATION_CALCULATION_FAILED

        # Get automated parallelization settings
        parallelization = get_pw_parallelization_parameters(calculation, **self.ctx.automatic_parallelization)

        # Note: don't do this at home, we are losing provenance here. This should be done by a calculation function
        node = orm.Dict(dict=parallelization).store()
        self.out('automatic_parallelization', node)
        self.report('results of automatic parallelization in {}<{}>'.format(node.__class__.__name__, node.pk))

        options = self.ctx.inputs.metadata['options']
        base_resources = options.get('resources', {})
        goal_resources = parallelization['resources']

        scheduler = calculation.computer.get_scheduler()
        resources = create_scheduler_resources(scheduler, base_resources, goal_resources)

        cmdline = self.ctx.inputs.settings.get('cmdline', [])
        cmdline = cmdline_remove_npools(cmdline)
        cmdline.extend(['-nk', str(parallelization['npools'])])

        # Set the new cmdline setting and resource options
        self.ctx.inputs.settings['cmdline'] = cmdline
        self.ctx.inputs.metadata['options'] = update_mapping(options, {'resources': resources})

        # Remove the only initialization flag
        self.ctx.inputs.settings.pop('ONLY_INITIALIZATION')

        return

    def prepare_calculation(self):
        """Prepare the inputs for the next calculation.

        If a `restart_calc` has been set in the context, its `remote_folder` will be used as the `parent_folder` input
        for the next calculation and the `restart_mode` is set to `restart`. Otherwise, no `parent_folder` is used and
        `restart_mode` is set to `from_scratch`.
        """
        if self.ctx.restart_calc:
            self.ctx.inputs.parameters['CONTROL']['restart_mode'] = 'restart'
            self.ctx.inputs.parent_folder = self.ctx.restart_calc.outputs.remote_folder
        else:
            self.ctx.inputs.parameters['CONTROL']['restart_mode'] = 'from_scratch'
            self.ctx.inputs.pop('parent_folder', None)

    def _handle_calculation_sanity_checks(self, calculation):
        """Perform sanity checks on the current `calculation` which has finished successfully according to the parser.

        Verify that the occupation of the last band is below a certain threshold, unless `occupations` was explicitly
        set to `fixed` in the input parameters. If this is violated, the calculation used too few bands and cannot be
        trusted. The number of bands is increased and the calculation is restarted, starting from the last.
        """
        from aiida_quantumespresso.utils.bands import get_highest_occupied_band

        # Only skip the check on the highest band occupation if `occupations` was explicitly set to `fixed`.
        if calculation.inputs.parameters.get_attribute('SYSTEM', {}).get('occupations', None) == 'fixed':
            return

        try:
            bands = calculation.outputs.output_band
            get_highest_occupied_band(bands)
        except ValueError as exception:
            self.report('calculation<{}> run with smearing and highest band is occupied'.format(calculation.pk))
            self.report('BandsData<{}> has invalid occupations: {}'.format(bands.pk, exception))
            return self._handle_insufficient_bands(calculation)

    def report_error_handled(self, calculation, action):
        """Report an action taken for a calculation that has failed.

        This should be called in a registered error handler if its condition is met and an action was taken.

        :param calculation: the failed calculation node
        :param action: a string message with the action taken
        """
        arguments = [calculation.process_label, calculation.pk, calculation.exit_status, calculation.exit_message]
        self.report('{}<{}> failed with exit status {}: {}'.format(*arguments))
        self.report('Action taken: {}'.format(action))


@register_error_handler(PwBaseWorkChain, 600)
def _handle_unrecoverable_failure(self, calculation):
    """Handle calculations with an exit status below 400 which are unrecoverable, so abort the work chain."""
    if calculation.exit_status < 400:
        self.report_error_handled(calculation, 'unrecoverable error, aborting...')
        return ErrorHandlerReport(True, True, self.exit_codes.ERROR_UNRECOVERABLE_FAILURE)


@register_error_handler(PwBaseWorkChain, 580)
def _handle_out_of_walltime(self, calculation):
    """Handle `ERROR_OUT_OF_WALLTIME` exit code: calculation shut down neatly and we can simply restart."""
    if calculation.exit_status == PwCalculation.spec().exit_codes.ERROR_OUT_OF_WALLTIME.status:
        try:
            self.ctx.inputs.structure = calculation.outputs.output_structure
        except exceptions.NotExistent:
            self.ctx.restart_calc = calculation
            self.report_error_handled(calculation, 'simply restart from the last calculation')
        else:
            self.ctx.restart_calc = None
            self.report_error_handled(calculation, 'out of walltime: structure changed so restarting from scratch')

        return ErrorHandlerReport(True, True)


@register_error_handler(PwBaseWorkChain, 570)
def _handle_vcrelax_converged_except_final_scf(self, calculation):
    """Handle `ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF` exit code.

    Convergence reached in `vc-relax` except thresholds exceeded in final scf: consider as converged.
    """
    exit_code_labels = [
        'ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF',
    ]

    if calculation.exit_status in PwCalculation.get_exit_statuses(exit_code_labels):
        self.ctx.is_finished = True
        self.ctx.restart_calc = calculation
        action = 'ionic convergence thresholds met except in final scf: consider structure relaxed.'
        self.report_error_handled(calculation, action)
        self.results()  # Call the results method to attach the output nodes
        return ErrorHandlerReport(True, True, self.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF)


@register_error_handler(PwBaseWorkChain, 560)
def _handle_relax_recoverable_ionic_convergence_error(self, calculation):
    """Handle various exit codes for recoverable `vc-relax` or `relax` calculations with failed ionic convergence.

    These exit codes signify that the ionic convergence thresholds were not met, but the output structure is usable, so
    the solution is to simply restart from scratch but from the output structure.
    """
    exit_code_labels = [
        'ERROR_IONIC_CONVERGENCE_NOT_REACHED',
        'ERROR_IONIC_CYCLE_EXCEEDED_NSTEP',
        'ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE',
        'ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE',
    ]

    if calculation.exit_status in PwCalculation.get_exit_statuses(exit_code_labels):
        self.ctx.restart_calc = None
        self.ctx.inputs.structure = calculation.outputs.output_structure
        action = 'no ionic convergence but clean shutdown: restarting from scratch but using output structure.'
        self.report_error_handled(calculation, action)
        return ErrorHandlerReport(True, True)


@register_error_handler(PwBaseWorkChain, 550)
def _handle_relax_recoverable_electronic_convergence_error(self, calculation):
    """Handle various exit codes for recoverable `vc-relax` or `relax` calculations with failed electronic convergence.

    These exit codes signify that the electronic convergence thresholds were not met, but the output structure is
    usable, so the solution is to simply restart from scratch but from the output structure.
    """
    exit_code_labels = [
        'ERROR_IONIC_CYCLE_ELECTRONIC_CONVERGENCE_NOT_REACHED',
        'ERROR_IONIC_CONVERGENCE_REACHED_FINAL_SCF_FAILED',
    ]

    if calculation.exit_status in PwCalculation.get_exit_statuses(exit_code_labels):

        factor = self.defaults.delta_factor_mixing_beta
        mixing_beta = self.ctx.inputs.parameters.get('ELECTRONS', {}).get('mixing_beta', self.defaults.qe.mixing_beta)
        mixing_beta_new = mixing_beta * factor

        self.ctx.restart_calc = None
        self.ctx.inputs.parameters.setdefault('ELECTRONS', {})['mixing_beta'] = mixing_beta_new
        self.ctx.inputs.structure = calculation.outputs.output_structure
        action = 'no electronic convergence but clean shutdown: reduced beta mixing from {} to {} restarting from ' \
                 'scratch but using output structure.'.format(mixing_beta, mixing_beta_new)
        self.report_error_handled(calculation, action)
        return ErrorHandlerReport(True, True)


@register_error_handler(PwBaseWorkChain, 410)
def _handle_electronic_convergence_not_achieved(self, calculation):
    """Handle `ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED`: decrease the mixing beta and restart from scratch."""
    if calculation.exit_status == PwCalculation.spec().exit_codes.ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED.status:
        factor = self.defaults.delta_factor_mixing_beta
        mixing_beta = self.ctx.inputs.parameters.get('ELECTRONS', {}).get('mixing_beta', self.defaults.qe.mixing_beta)
        mixing_beta_new = mixing_beta * factor

        self.ctx.restart_calc = calculation
        self.ctx.inputs.parameters.setdefault('ELECTRONS', {})['mixing_beta'] = mixing_beta_new

        action = 'reduced beta mixing from {} to {} and restarting from the last ' \
            'calculation'.format(mixing_beta, mixing_beta_new)
        self.report_error_handled(calculation, action)
        return ErrorHandlerReport(True, True)


@register_error_handler(PwBaseWorkChain)
def _handle_insufficient_bands(self, calculation):
    """Handle successfully converged calculation with too few bands, so increase them and restart from scratch."""
    nbnd_cur = calculation.outputs.output_parameters.get_dict()['number_of_bands']
    nbnd_new = nbnd_cur + max(int(nbnd_cur * self.defaults.delta_factor_nbnd), self.defaults.delta_minimum_nbnd)

    self.ctx.inputs.parameters.setdefault('SYSTEM', {})['nbnd'] = nbnd_new
    self.ctx.restart_calc = None

    self.report('{}<{}> had insufficient bands'.format(calculation.process_label, calculation.pk))
    self.report('Action taken: increased the number of bands to {} and restarting from scratch'.format(nbnd_new))
    return ErrorHandlerReport(True, True)
