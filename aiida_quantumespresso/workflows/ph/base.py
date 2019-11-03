# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO ph.x calculation with automated error handling and restarts."""
from __future__ import absolute_import

from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import while_
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.common.workchain.utils import ErrorHandlerReport, register_error_handler
from aiida_quantumespresso.common.workchain.base.restart import BaseRestartWorkChain

PhCalculation = CalculationFactory('quantumespresso.ph')
PwCalculation = CalculationFactory('quantumespresso.pw')


class PhBaseWorkChain(BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO ph.x calculation with automated error handling and restarts."""

    _calculation_class = PhCalculation

    defaults = AttributeDict({
        'delta_factor_max_seconds': 0.95,
        'delta_factor_alpha_mix': 0.90,
        'alpha_mix': 0.70,
    })

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super(PhBaseWorkChain, cls).define(spec)
        spec.expose_inputs(PhCalculation, namespace='ph')
        spec.input('only_initialization', valid_type=orm.Bool, default=orm.Bool(False))
        spec.outline(
            cls.setup,
            cls.validate_parameters,
            cls.validate_resources,
            while_(cls.should_run_calculation)(
                cls.prepare_calculation,
                cls.run_calculation,
                cls.inspect_calculation,
            ),
            cls.results,
        )
        spec.expose_outputs(PwCalculation, exclude=('retrieved_folder',))
        spec.exit_code(204, 'ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED',
            message='The `metadata.options` did not specify both `resources.num_machines` and `max_wallclock_seconds`.')
        spec.exit_code(300, 'ERROR_UNRECOVERABLE_FAILURE',
            message='The calculation failed with an unrecoverable error.')

    def setup(self):
        """Call the `setup` of the `BaseRestartWorkChain` and then create the inputs dictionary in `self.ctx.inputs`.

        This `self.ctx.inputs` dictionary will be used by the `BaseRestartWorkChain` to submit the calculations in the
        internal loop.
        """
        super(PhBaseWorkChain, self).setup()
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

    def prepare_calculation(self):
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
        self.report('Action taken: {}'.format(action))


@register_error_handler(PhBaseWorkChain, 600)
def _handle_unrecoverable_failure(self, calculation):
    """Jandle calculations with an exit status below 400 which are unrecoverable, so abort the work chain."""
    if calculation.exit_status < 400:
        self.report_error_handled(calculation, 'unrecoverable error, aborting...')
        return ErrorHandlerReport(True, True, self.exit_codes.ERROR_UNRECOVERABLE_FAILURE)


@register_error_handler(PhBaseWorkChain, 580)
def _handle_out_of_walltime(self, calculation):
    """Handle `ERROR_OUT_OF_WALLTIME` exit code: calculation shut down neatly and we can simply restart."""
    if calculation.exit_status == PhCalculation.exit_codes.ERROR_OUT_OF_WALLTIME.status:
        self.ctx.restart_calc = calculation
        self.report_error_handled(calculation, 'simply restart from the last calculation')
        return ErrorHandlerReport(True, True)


@register_error_handler(PhBaseWorkChain, 410)
def _handle_convergence_not_achieved(self, calculation):
    """Handle `ERROR_CONVERGENCE_NOT_REACHED` exit code: decrease the mixing beta and restart from scratch."""
    if calculation.exit_status == PwCalculation.exit_codes.ERROR_CONVERGENCE_NOT_REACHED.status:
        factor = self.defaults.delta_factor_alpha_mix
        alpha_mix = self.ctx.inputs.parameters.get('INPUTPH', {}).get('alpha_mix(1)', self.defaults.alpha_mix)
        alpha_mix_new = alpha_mix * factor

        self.ctx.restart_calc = calculation
        self.ctx.inputs.parameters.setdefault('INPUTPH', {})['alpha_mix(1)'] = alpha_mix_new

        action = 'reduced alpha_mix from {} to {} and restarting'.format(alpha_mix, alpha_mix_new)
        self.report_error_handled(calculation, action)
        return ErrorHandlerReport(True, True)
