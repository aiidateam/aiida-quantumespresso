"""Workchain to run a Quantum ESPRESSO q2r.x calculation with automated error handling and restarts."""

from aiida.common import AttributeDict
from aiida.engine import BaseRestartWorkChain, while_
from aiida.plugins import CalculationFactory

Q2rCalculation = CalculationFactory('quantumespresso.q2r')


class Q2rBaseWorkChain(BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO q2r.x calculation with automated error handling and restarts."""

    _process_class = Q2rCalculation

    @classmethod
    def define(cls, spec):
        """Define the process specification."""

        super().define(spec)
        spec.expose_inputs(Q2rCalculation, namespace='q2r')
        spec.expose_outputs(Q2rCalculation)
        spec.outline(
            cls.setup,
            while_(cls.should_run_process)(
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )
        spec.exit_code(
            300,
            'ERROR_UNRECOVERABLE_FAILURE',
            message='[deprecated] The calculation failed with an unrecoverable error.',
        )

    def setup(self):
        """Call the `setup` of the `BaseRestartWorkChain` and then create the inputs dictionary in `self.ctx.inputs`.

        This `self.ctx.inputs` dictionary will be used by the `BaseRestartWorkChain` to submit the calculations in the
        internal loop.
        """
        super().setup()
        self.ctx.restart_calc = None
        self.ctx.inputs = AttributeDict(self.exposed_inputs(Q2rCalculation, 'q2r'))

    def report_error_handled(self, calculation, action):
        """Report an action taken for a calculation that has failed.

        This should be called in a registered error handler if its condition is met and an action was taken.

        :param calculation: the failed calculation node
        :param action: a string message with the action taken
        """
        arguments = [
            calculation.process_label,
            calculation.pk,
            calculation.exit_status,
            calculation.exit_message,
        ]
        self.report('{}<{}> failed with exit status {}: {}'.format(*arguments))
        self.report(f'Action taken: {action}')
