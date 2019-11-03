# -*- coding: utf-8 -*-
"""Base implementation of `WorkChain` class that implements a simple automated restart mechanism for calculations."""
from __future__ import absolute_import

from aiida import orm
from aiida.common import exceptions
from aiida.common.lang import override
from aiida.engine import CalcJob, WorkChain, ExitCode, ToContext, append_
from aiida.plugins.entry_point import get_entry_point_names, load_entry_point

from aiida_quantumespresso.common.exceptions import UnexpectedCalculationFailure
from aiida_quantumespresso.common.workchain.utils import ErrorHandlerReport


class BaseRestartWorkChain(WorkChain):
    """Base restart work chain.

    This work chain serves as the starting point for more complex work chains that will be designed to run a calculation
    that might need multiple restarts to come to a successful end. These restarts may be necessary because a single
    calculation run is not sufficient to achieve a fully converged result, or certain errors maybe encountered which
    are recoverable.

    This work chain implements the most basic functionality to achieve this goal. It will launch calculations,
    restarting until it is completed successfully or the maximum number of iterations is reached. It can recover from
    errors through error handlers that can be attached dynamically through the `register_error_handler` decorator.

    The idea is to sub class this work chain and leverage the generic error handling that is implemented in the few
    outline methods. The minimally required outline would look something like the following::

        cls.setup
        while_(cls.should_run_calculation)(
            cls.run_calculation,
            cls.inspect_calculation,
        )

    Each of these methods can of course be overriden but they should be general enough to fit most calculation cycles.
    The `run_calculation` method will take the inputs for the calculation process from the context under the key
    `inputs`. The user should therefore make sure that before the `run_calculation` method is called, that the to be
    used inputs are stored under `self.ctx.inputs`. One can update the inputs based on the results from a prior
    calculation by calling an outline method just before the `run_calculation` step, for example::

        cls.setup
        while_(cls.should_run_calculation)(
            cls.prepare_calculation,
            cls.run_calculation,
            cls.inspect_calculation,
        )

    Where in the `prepare_calculation` method, the inputs dictionary at `self.ctx.inputs` is updated before the next
    calculation will be run with those inputs.

    The `_calculation_class` attribute should be set to the `CalcJob` class that should be run in the loop.
    """

    _verbose = False
    _calculation_class = None
    _error_handler_entry_point = None

    def __init__(self, *args, **kwargs):
        """Construct the instance."""
        super(BaseRestartWorkChain, self).__init__(*args, **kwargs)

        if self._calculation_class is None or not issubclass(self._calculation_class, CalcJob):
            raise ValueError('no valid CalcJob class defined for `_calculation_class` attribute')

        self._load_error_handlers()

    @override
    def load_instance_state(self, saved_state, load_context):
        """Load the process instance from a saved state.

        :param saved_state: saved state of existing process instance
        :param load_context: context for loading instance state
        """
        super(BaseRestartWorkChain, self).load_instance_state(saved_state, load_context)
        self._load_error_handlers()

    def _load_error_handlers(self):
        """Load the error handlers defined through entry points, if any."""
        if self._error_handler_entry_point is not None:
            for entry_point_name in get_entry_point_names(self._error_handler_entry_point):
                try:
                    load_entry_point(self._error_handler_entry_point, entry_point_name)
                    self.logger.info(
                        "loaded the '%s' entry point for the '%s' error handlers category", entry_point_name,
                        self._error_handler_entry_point
                    )
                except exceptions.EntryPointError as exception:
                    self.logger.warning(
                        "failed to load the '%s' entry point for the '%s' error handlers: %s", entry_point_name,
                        self._error_handler_entry_point, exception
                    )

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super(BaseRestartWorkChain, cls).define(spec)
        spec.input('max_iterations', valid_type=orm.Int, default=orm.Int(5),
            help='Maximum number of iterations the work chain will restart the calculation to finish successfully.')
        spec.input('clean_workdir', valid_type=orm.Bool, default=orm.Bool(False),
            help='If `True`, work directories of all called calculation will be cleaned at the end of execution.')
        spec.exit_code(101, 'ERROR_MAXIMUM_ITERATIONS_EXCEEDED',
            message='The maximum number of iterations was exceeded.')
        spec.exit_code(102, 'ERROR_SECOND_CONSECUTIVE_UNHANDLED_FAILURE',
            message='The calculation failed for an unknown reason, twice in a row.')

    def setup(self):
        """Initialize context variables that are used during the logical flow of the `BaseRestartWorkChain`."""
        self.ctx.calc_name = self._calculation_class.__name__
        self.ctx.unexpected_failure = False
        self.ctx.restart_calc = None
        self.ctx.is_finished = False
        self.ctx.iteration = 0

    def should_run_calculation(self):
        """Return whether a new calculation should be run.

        This is the case as long as the last calculation has not finished successfully and the maximum number of
        restarts has not yet been exceeded.
        """
        return not self.ctx.is_finished and self.ctx.iteration < self.inputs.max_iterations.value

    def run_calculation(self):
        """Run the next calculation, taking the input dictionary from the context at `self.ctx.inputs`."""
        from aiida_quantumespresso.utils.mapping import prepare_process_inputs

        self.ctx.iteration += 1

        try:
            unwrapped_inputs = self.ctx.inputs
        except AttributeError:
            raise AttributeError('no calculation input dictionary was defined in `self.ctx.inputs`')

        # Set the `CALL` link label
        unwrapped_inputs['metadata']['call_link_label'] = 'iteration_{:02d}'.format(self.ctx.iteration)

        inputs = prepare_process_inputs(self._calculation_class, unwrapped_inputs)
        calculation = self.submit(self._calculation_class, **inputs)

        # Add a new empty list to the `errors_handled` extra. If any errors handled registered through the
        # `register_error_handler` decorator return an `ErrorHandlerReport`, their name will be appended to that list.
        errors_handled = self.node.get_extra('errors_handled', [])
        errors_handled.append([])
        self.node.set_extra('errors_handled', errors_handled)

        self.report('launching {}<{}> iteration #{}'.format(self.ctx.calc_name, calculation.pk, self.ctx.iteration))

        return ToContext(calculations=append_(calculation))

    def inspect_calculation(self):
        """Analyse the results of the previous calculation and call the error handlers when necessary."""
        calculation = self.ctx.calculations[self.ctx.iteration - 1]

        # Done: successful completion of last calculation
        if calculation.is_finished_ok:

            # Perform an optional sanity check. If it returns an `ExitCode` this means an unrecoverable situation was
            # detected and the work chain should be aborted. If it returns `False`, the sanity check detected a problem
            # but has handled the problem and we should restart the cycle.
            handler = self._handle_calculation_sanity_checks(calculation)  # pylint: disable=assignment-from-no-return

            if (isinstance(handler, ErrorHandlerReport) and
                    handler.exit_code is not None and handler.exit_code.status != 0):
                # Sanity check returned a handler with an exit code that is non-zero, so we abort
                self.report('{}<{}> finished successfully, but sanity check detected unrecoverable problem'.format(
                    self.ctx.calc_name, calculation.pk))
                return handler.exit_code

            if isinstance(handler, ErrorHandlerReport):
                # Reset the `unexpected_failure` since we are restarting the calculation loop
                self.ctx.unexpected_failure = False
                self.report('{}<{}> finished successfully, but sanity check failed, restarting'.format(
                    self.ctx.calc_name, calculation.pk))
                return

            self.report('{}<{}> completed successfully'.format(self.ctx.calc_name, calculation.pk))
            self.ctx.restart_calc = calculation
            self.ctx.is_finished = True
            return

        # Unexpected: calculation was killed or an exception occurred, trigger unexpected failure handling
        if calculation.is_excepted or calculation.is_killed:
            return self._handle_unexpected_failure(calculation)

        # Failed: here the calculation is `Finished` but has a non-zero exit status, initiate the error handling
        try:
            exit_code = self._handle_calculation_failure(calculation)
        except UnexpectedCalculationFailure as exception:
            exit_code = self._handle_unexpected_failure(calculation, exception)

        # If the exit code returned actually has status `0` that means we consider the calculation as successful
        if isinstance(exit_code, ExitCode) and exit_code.status == 0:
            self.ctx.is_finished = True

        return exit_code

    def results(self):
        """Attach the outputs specified in the output specification from the last completed calculation."""
        calculation = self.ctx.calculations[self.ctx.iteration - 1]

        # We check the `is_finished` attribute of the work chain and not the successfulness of the last calculation
        # because the error handlers in the last iteration can have qualified a "failed" calculation as satisfactory
        # for the outcome of the work chain and so have marked it as `is_finished=True`.
        if not self.ctx.is_finished and self.ctx.iteration >= self.inputs.max_iterations.value:
            # Abort: exceeded maximum number of retries
            self.report('reached the maximum number of iterations {}: last ran {}<{}>'.format(
                self.inputs.max_iterations.value, self.ctx.calc_name, calculation.pk))
            return self.exit_codes.ERROR_MAXIMUM_ITERATIONS_EXCEEDED

        self.report('work chain completed after {} iterations'.format(self.ctx.iteration))

        for name, port in self.spec().outputs.items():

            try:
                node = calculation.get_outgoing(link_label_filter=name).one().node
            except ValueError:
                if port.required:
                    self.report("required output '{}' was not an output of {}<{}>".format(
                        name, self.ctx.calc_name, calculation.pk))
            else:
                self.out(name, node)
                if self._verbose:
                    self.report("attaching the node {}<{}> as '{}'".format(node.__class__.__name__, node.pk, name))

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super(BaseRestartWorkChain, self).on_terminated()

        if self.inputs.clean_workdir.value is False:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = []

        for called_descendant in self.node.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                try:
                    called_descendant.outputs.remote_folder._clean()  # pylint: disable=protected-access
                    cleaned_calcs.append(str(called_descendant.pk))
                except (IOError, OSError, KeyError):
                    pass

        if cleaned_calcs:
            self.report('cleaned remote folders of calculations: {}'.format(' '.join(cleaned_calcs)))

    def _handle_calculation_sanity_checks(self, calculation):
        """Perform a sanity check of a calculation that finished ok.

        Calculations that were marked as successful by the parser may still have produced outputs that do not make sense
        but were not detected by the code and so were not highlighted as warnings or errors. The consistency of the
        outputs can be checked here. If an unrecoverable problem is found, the function should return the appropriate
        exit code to abort the work chain. If the probem can be fixed with a restart calculation, this function should
        adapt the inputs as an error handler would and return `False`. This will signal to the work chain that a new
        calculation should be started. If `None` is returned, the work chain assumes that the outputs produced by the
        calculation are good and nothing will be done.

        :param calculation: the calculation whose outputs should be checked for consistency
        :return: `ErrorHandlerReport` if a new calculation should be launched or abort if it includes an exit code
        """

    def _handle_calculation_failure(self, calculation):
        """Call the attached error handlers if any to attempt to correct the cause of the calculation failure.

        The registered error handlers will be called in order based on their priority until a handler returns a report
        that instructs to break. If the last executed error handler defines an exit code, that will be returned to
        instruct the work chain to abort. Otherwise the work chain will continue the cycle.

        :param calculation: the calculation that finished with a non-zero exit status
        :return: `ExitCode` if the work chain is to be aborted
        :raises `UnexpectedCalculationFailure`: if no error handlers were registered or no errors were handled.
        """
        is_handled = False
        handler_report = None

        if not hasattr(self, '_error_handlers') or not self._error_handlers:
            raise UnexpectedCalculationFailure('no calculation error handlers were registered')

        # Sort the handlers with a priority defined, based on their priority in reverse order
        handlers = [handler for handler in self._error_handlers if handler.priority]
        handlers = sorted(handlers, key=lambda x: x.priority, reverse=True)

        for handler in handlers:

            handler_report = handler.method(self, calculation)

            # If at least one error is handled, we consider the calculation failure handled.
            if handler_report and handler_report.is_handled:
                self.ctx.unexpected_failure = False
                is_handled = True

            # After certain error handlers, we may want to skip all other error handling
            if handler_report and handler_report.do_break:
                break

        # If none of the executed error handlers reported that they handled an error, the failure reason is unknown
        if not is_handled:
            raise UnexpectedCalculationFailure('calculation failure was not handled')

        # The last called error handler may not necessarily have returned a handler report
        if handler_report:
            return handler_report.exit_code

        return

    def _handle_unexpected_failure(self, calculation, exception=None):
        """Handle an unexpected failure.

        This occurs when a calculation excepted, was killed or finished with a non-zero exit status but no errors were
        handled. If this is the second consecutive unexpected failure the work chain is aborted.

        :param calculation: the calculation that failed in an unexpected way
        :param exception: optional exception or error message to log to the report
        :return: `ExitCode` if this is the second consecutive unexpected failure
        """
        if exception:
            self.report('{}'.format(exception))

        if self.ctx.unexpected_failure:
            self.report('failure of {}<{}> could not be handled for the second consecutive time'.format(
                self.ctx.calc_name, calculation.pk))
            return self.exit_codes.ERROR_SECOND_CONSECUTIVE_UNHANDLED_FAILURE

        self.ctx.unexpected_failure = True
        self.report('failure of {}<{}> could not be handled, restarting once more'.format(
            self.ctx.calc_name, calculation.pk))
