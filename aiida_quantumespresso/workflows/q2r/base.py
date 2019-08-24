# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO q2r.x calculation with automated error handling and restarts."""
from __future__ import absolute_import

from aiida.common import AttributeDict
from aiida.engine import while_
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.common.workchain.base.restart import BaseRestartWorkChain
from aiida_quantumespresso.common.workchain.utils import register_error_handler, ErrorHandlerReport

Q2rCalculation = CalculationFactory('quantumespresso.q2r')


class Q2rBaseWorkChain(BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO q2r.x calculation with automated error handling and restarts."""

    _calculation_class = Q2rCalculation

    @classmethod
    def define(cls, spec):
        # yapf: disable
        super(Q2rBaseWorkChain, cls).define(spec)
        spec.expose_inputs(Q2rCalculation, namespace='q2r')
        spec.expose_outputs(Q2rCalculation)
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
        super(Q2rBaseWorkChain, self).setup()
        self.ctx.inputs = AttributeDict(self.exposed_inputs(Q2rCalculation, 'q2r'))

    def report_error_handled(self, calculation, action):
        """Report an action taken for a calculation that has failed.

        This should be called in a registered error handler if its condition is met and an action was taken.

        :param calculation: the failed calculation node
        :param action: a string message with the action taken
        """
        arguments = [calculation.process_label, calculation.pk, calculation.exit_status, calculation.exit_message]
        self.report('{}<{}> failed with exit status {}: {}'.format(*arguments))
        self.report('Action taken: {}'.format(action))


@register_error_handler(Q2rBaseWorkChain, 600)
def _handle_unrecoverable_failure(self, calculation):
    """Calculations with an exit status below 400 are unrecoverable, so abort the work chain."""
    if calculation.exit_status < 400:
        self.report_error_handled(calculation, 'unrecoverable error, aborting...')
        return ErrorHandlerReport(True, True, self.exit_codes.ERROR_UNRECOVERABLE_FAILURE)
