# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO ph.x calculation with automated error handling and restarts."""
from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import while_, BaseRestartWorkChain, process_handler, ProcessHandlerReport
from aiida.plugins import CalculationFactory

PhCalculation = CalculationFactory('quantumespresso.ph')
PwCalculation = CalculationFactory('quantumespresso.pw')


class PhBaseWorkChain(BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO ph.x calculation with automated error handling and restarts."""

    _process_class = PhCalculation

    defaults = AttributeDict({
        'delta_factor_max_seconds': 0.95,
        'delta_factor_alpha_mix': 0.90,
        'alpha_mix': 0.70,
    })

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.expose_inputs(PhCalculation, namespace='ph')
        spec.input('only_initialization', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.outline(
            cls.setup,
            cls.validate_parameters,
            cls.validate_resources,
            while_(cls.should_run_process)(
                cls.prepare_process,
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )
        spec.expose_outputs(PwCalculation, exclude=('retrieved_folder',))
        spec.exit_code(204, 'ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED',
            message='The `metadata.options` did not specify both `resources.num_machines` and `max_wallclock_seconds`.')
        spec.exit_code(300, 'ERROR_UNRECOVERABLE_FAILURE',
            message='The calculation failed with an unrecoverable error.')
        # yapf: enable

    def setup(self):
        """Call the `setup` of the `BaseRestartWorkChain` and then create the inputs dictionary in `self.ctx.inputs`.

        This `self.ctx.inputs` dictionary will be used by the `BaseRestartWorkChain` to submit the calculations in the
        internal loop.
        """
        super().setup()
        self.ctx.restart_calc = None
        self.ctx.inputs = AttributeDict(self.exposed_inputs(PhCalculation, 'ph'))

    def validate_parameters(self):
        """Validate inputs that might depend on each other and cannot be validated by the spec."""
        self.ctx.inputs.parameters = self.ctx.inputs.parameters.get_dict()
        self.ctx.inputs.settings = self.ctx.inputs.settings.get_dict() if 'settings' in self.ctx.inputs else {}

        self.ctx.inputs.parameters.setdefault('INPUTPH', {})
        self.ctx.inputs.parameters['INPUTPH']['recover'] = 'parent_folder' in self.ctx.inputs

        if self.inputs.only_initialization.value:
            self.ctx.inputs.settings['ONLY_INITIALIZATION'] = True

    def validate_resources(self):
        """Validate the inputs related to the resources.

        The `metadata.options` should at least contain the options `resources` and `max_wallclock_seconds`, where
        `resources` should define the `num_machines`.
        """
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
        self.ctx.inputs.parameters['INPUTPH']['max_seconds'] = max_seconds

    def prepare_process(self):
        """Prepare the inputs for the next calculation.

        If a `restart_calc` has been set in the context, its `remote_folder` will be used as the `parent_folder` input
        for the next calculation and the `restart_mode` is set to `restart`.
        """
        if self.ctx.restart_calc:
            self.ctx.inputs.parameters['INPUTPH']['recover'] = True
            self.ctx.inputs.parent_folder = self.ctx.restart_calc.outputs.remote_folder

    def report_error_handled(self, calculation, action):
        """Report an action taken for a calculation that has failed.

        This should be called in a registered error handler if its condition is met and an action was taken.

        :param calculation: the failed calculation node
        :param action: a string message with the action taken
        """
        arguments = [calculation.process_label, calculation.pk, calculation.exit_status, calculation.exit_message]
        self.report('{}<{}> failed with exit status {}: {}'.format(*arguments))
        self.report(f'Action taken: {action}')

    @process_handler(priority=600)
    def handle_unrecoverable_failure(self, node):
        """Handle calculations with an exit status below 400 which are unrecoverable, so abort the work chain."""
        if node.is_failed and node.exit_status < 400:
            self.report_error_handled(node, 'unrecoverable error, aborting...')
            return ProcessHandlerReport(True, self.exit_codes.ERROR_UNRECOVERABLE_FAILURE)

    @process_handler(priority=580, exit_codes=PhCalculation.exit_codes.ERROR_OUT_OF_WALLTIME)
    def handle_out_of_walltime(self, node):
        """Handle `ERROR_OUT_OF_WALLTIME` exit code: calculation shut down neatly and we can simply restart."""
        self.ctx.restart_calc = node
        self.report_error_handled(node, 'simply restart from the last calculation')
        return ProcessHandlerReport(True)

    @process_handler(priority=410, exit_codes=PhCalculation.exit_codes.ERROR_CONVERGENCE_NOT_REACHED)
    def handle_convergence_not_achieved(self, node):
        """Handle `ERROR_CONVERGENCE_NOT_REACHED` exit code: decrease the mixing beta and restart from scratch."""
        factor = self.defaults.delta_factor_alpha_mix
        alpha_mix = self.ctx.inputs.parameters.get('INPUTPH', {}).get('alpha_mix(1)', self.defaults.alpha_mix)
        alpha_mix_new = alpha_mix * factor

        self.ctx.restart_calc = node
        self.ctx.inputs.parameters.setdefault('INPUTPH', {})['alpha_mix(1)'] = alpha_mix_new

        action = f'reduced alpha_mix from {alpha_mix} to {alpha_mix_new} and restarting'
        self.report_error_handled(node, action)
        return ProcessHandlerReport(True)
