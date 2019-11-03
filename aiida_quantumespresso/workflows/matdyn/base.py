# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO matdyn.x calculation with automated error handling and restarts."""
from __future__ import absolute_import

from aiida.common import AttributeDict
from aiida.engine import while_
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.common.workchain.base.restart import BaseRestartWorkChain
from aiida_quantumespresso.common.workchain.utils import register_error_handler, ErrorHandlerReport

MatdynCalculation = CalculationFactory('quantumespresso.matdyn')


class MatdynBaseWorkChain(BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO matdyn.x calculation with automated error handling and restarts."""

    _calculation_class = MatdynCalculation

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super(MatdynBaseWorkChain, cls).define(spec)
        spec.expose_inputs(MatdynCalculation, namespace='matdyn')
        spec.expose_outputs(MatdynCalculation)
        spec.outline(
            cls.setup,
            while_(cls.should_run_calculation)(
                cls.run_calculation,
                cls.inspect_calculation,
            ),
            cls.results,
        )
        spec.exit_code(300, 'ERROR_UNRECOVERABLE_FAILURE',
            message='The calculation failed with an unrecoverable error.')

    def setup(self):
        """Call the `setup` of the `BaseRestartWorkChain` and then create the inputs dictionary in `self.ctx.inputs`.

        This `self.ctx.inputs` dictionary will be used by the `BaseRestartWorkChain` to submit the calculations in the
        internal loop.
        """
        super(MatdynBaseWorkChain, self).setup()
        self.ctx.inputs = AttributeDict(self.exposed_inputs(MatdynCalculation, 'matdyn'))

    def report_error_handled(self, calculation, action):
        """Report an action taken for a calculation that has failed.

        This should be called in a registered error handler if its condition is met and an action was taken.

        :param calculation: the failed calculation node
        :param action: a string message with the action taken
        """
        arguments = [calculation.process_label, calculation.pk, calculation.exit_status, calculation.exit_message]
        self.report('{}<{}> failed with exit status {}: {}'.format(*arguments))
        self.report('Action taken: {}'.format(action))


@register_error_handler(MatdynBaseWorkChain, 600)
def _handle_unrecoverable_failure(self, calculation):
    """Handle calculations with an exit status below 400 which are unrecoverable, so abort the work chain."""
    if calculation.exit_status < 400:
        self.report_error_handled(calculation, 'unrecoverable error, aborting...')
        return ErrorHandlerReport(True, True, self.exit_codes.ERROR_UNRECOVERABLE_FAILURE)
