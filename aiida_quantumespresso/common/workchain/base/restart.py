# -*- coding: utf-8 -*-
from aiida.common.datastructures import calc_states
from aiida.common.exceptions import LoadingPluginFailed, MissingPluginError
from aiida.orm.calculation import JobCalculation
from aiida.orm.data.base import Bool, Int
from aiida.work.workchain import WorkChain, ToContext, append_
from aiida_quantumespresso.common.exceptions import UnexpectedCalculationFailure
from aiida_quantumespresso.common.pluginloader import get_plugin, get_plugins


class BaseRestartWorkChain(WorkChain):
    """
    Base restart workchain

    This workchain serves as the starting point for more complex workchains that will be designed to
    run a calculation that might need multiple restarts to come to a successful end. These restarts
    may be necessary because a single calculation run is not sufficient to achieve a fully converged
    result, or certain errors maybe encountered which are recoverable.

    This workchain implements the most basic functionality to achieve this goal. It will launch
    calculations, restarting until it is completed successfully or the maximum number of iterations
    is reached. It can recover from simple submission failures and will handle any other errors through
    any of the error handlers that have been provided.

    The idea is to subclass this workchain and leverage the generic error handling that is implemented
    in the few outline methods. Part of the suggested outline would look something like the following::

        cls.setup
        while_(cls.should_run_calculation)(
            cls.run_calculation,
            cls.inspect_calculation,
        )

    Each of these methods can of course be overriden but they should be general enough to fit most
    calculation cycles. The run_calculation method will take the inputs for the calculation process
    from the context under the key 'inputs'. The user should therefore make sure that before the
    run_calculation method is called, that the to be used inputs are stored under self.ctx.inputs.
    One can update the inputs based on the results from a prior calculation by calling an outline
    method just before the run_calculation step, for example::

        cls.setup
        while_(cls.should_run_calculation)(
            cls.prepare_calculation,
            cls.run_calculation,
            cls.inspect_calculation,
        )

    Where in the prepare_calculation method, the inputs dictionary at self.ctx.inputs is updated
    before the next calculation will be run with those inputs.
    """
    _verbose = False
    _calculation_class = None
    _error_handler_entry_point = None
    _expected_calculation_states = [calc_states.FINISHED, calc_states.FAILED, calc_states.SUBMISSIONFAILED]

    def __init__(self, *args, **kwargs):
        super(BaseRestartWorkChain, self).__init__(*args, **kwargs)

        if self._calculation_class is None or not issubclass(self._calculation_class, JobCalculation):
            raise ValueError('no valid JobCalculation class defined for _calculation_class attribute')

        # If an error handler entry point is defined, load them. If the plugin cannot be loaded log it and pass
        if self._error_handler_entry_point is not None:
            for plugin in get_plugins(self._error_handler_entry_point):
                try:
                    get_plugin(self._error_handler_entry_point, plugin)
                    self.logger.info("loaded the '{}' entry point for the '{}' error handlers category"
                        .format(plugin, self._error_handler_entry_point, plugin))
                except (LoadingPluginFailed, MissingPluginError):
                    self.logger.warning("failed to load the '{}' entry point for the '{}' error handlers"
                        .format(plugin, self._error_handler_entry_point))

        return

    @classmethod
    def define(cls, spec):
        super(BaseRestartWorkChain, cls).define(spec)
        spec.input('max_iterations', valid_type=Int, default=Int(5),
            help='maximum number of iterations the workchain will restart the calculation to finish successfully')
        spec.input('clean_workdir', valid_type=Bool, default=Bool(False),
            help='if True, work directories of all called calculation will be cleaned at the end of execution')
        spec.exit_code(100, 'ERROR_ITERATION_RETURNED_NO_CALCULATION',
            message='the run_calculation step did not successfully add a calculation node to the context')
        spec.exit_code(101, 'ERROR_MAXIMUM_ITERATIONS_EXCEEDED',
            message='the maximum number of iterations was exceeded')
        spec.exit_code(102, 'ERROR_UNEXPECTED_CALCULATION_STATE',
            message='the calculation finished with an unexpected calculation state')
        spec.exit_code(103, 'ERROR_SECOND_CONSECUTIVE_SUBMISSION_FAILURE',
            message='the calculation failed to submit, twice in a row')
        spec.exit_code(104, 'ERROR_SECOND_CONSECUTIVE_UNHANDLED_FAILURE',
            message='the calculation failed for an unknown reason, twice in a row')

    def setup(self):
        """
        Initialize context variables that are used during the logical flow of the BaseRestartWorkChain
        """
        self.ctx.calc_name = self._calculation_class.__name__
        self.ctx.unexpected_failure = False
        self.ctx.submission_failure = False
        self.ctx.restart_calc = None
        self.ctx.is_finished = False
        self.ctx.iteration = 0

        return

    def should_run_calculation(self):
        """
        Return whether a new calculation should be run, which is the case as long as the last
        calculation has not finished successfully and the maximum number of restarts has not yet
        been exceeded
        """
        return not self.ctx.is_finished and self.ctx.iteration < self.inputs.max_iterations.value

    def run_calculation(self):
        """
        Run a new calculation, taking the input dictionary from the context at self.ctx.inputs
        """
        self.ctx.iteration += 1

        try:
            unwrapped_inputs = self.ctx.inputs
        except AttributeError:
            raise ValueError('no calculation input dictionary was defined in self.ctx.inputs')

        inputs = self._prepare_process_inputs(self._calculation_class, unwrapped_inputs)
        calculation = self.submit(self._calculation_class, **inputs)

        self.report('launching {}<{}> iteration #{}'.format(self.ctx.calc_name, calculation.pk, self.ctx.iteration))

        return ToContext(calculations=append_(calculation))

    def inspect_calculation(self):
        """
        Analyse the results of the previous calculation, checking whether it finished successfully
        or if not troubleshoot the cause and handle the errors, or abort if unrecoverable error was found
        """
        try:
            calculation = self.ctx.calculations[self.ctx.iteration - 1]
        except IndexError:
            self.report('iteration {} finished without returning a {}'.format(self.ctx.iteration, self.ctx.calc_name))
            return self.exit_codes.ERROR_ITERATION_RETURNED_NO_CALCULATION

        exit_code = None

        # Done: successful completion of last calculation
        if calculation.is_finished_ok:
            self.report('{}<{}> completed successfully'.format(self.ctx.calc_name, calculation.pk))
            self.ctx.restart_calc = calculation
            self.ctx.is_finished = True

        # Abort: exceeded maximum number of retries
        elif self.ctx.iteration >= self.inputs.max_iterations.value:
            self.report('reached the maximumm number of iterations {}: last ran {}<{}>'
                .format(self.inputs.max_iterations.value, self.ctx.calc_name, calculation.pk))
            exit_code = self.exit_codes.ERROR_MAXIMUM_ITERATIONS_EXCEEDED

        # Abort: unexpected state of last calculation
        elif calculation.get_state() not in self._expected_calculation_states:
            self.report('unexpected state ({}) of {}<{}>'
                .format(calculation.get_state(), self.ctx.calc_name, calculation.pk))
            exit_code = self.exit_codes.ERROR_UNEXPECTED_CALCULATION_STATE

        # Retry or abort: submission failed, try to restart or abort
        elif calculation.get_state() in [calc_states.SUBMISSIONFAILED]:
            exit_code = self._handle_submission_failure(calculation)
            self.ctx.submission_failure = True

        # Retry or abort: calculation finished or failed
        elif calculation.get_state() in [calc_states.FINISHED, calc_states.FAILED]:

            # Calculation was at least submitted successfully, so we reset the flag
            self.ctx.submission_failure = False

            # Check output for problems independent on calculation state and that do not trigger parser warnings
            exit_code = self._handle_calculation_sanity_checks(calculation)

            # Calculation failed, try to salvage it or handle any unexpected failures
            if calculation.get_state() in [calc_states.FAILED]:
                try:
                    exit_code = self._handle_calculation_failure(calculation)
                except UnexpectedCalculationFailure as exception:
                    exit_code = self._handle_unexpected_failure(calculation, exception)
                    self.ctx.unexpected_failure = True

            # Calculation finished: but did not finish ok, simply try to restart from this calculation
            else:
                self.ctx.unexpected_failure = False
                self.ctx.restart_calc = calculation
                self.report('calculation terminated without errors but did not complete successfully, restarting')

        return exit_code

    def results(self):
        """
        Attach the outputs specified in the output specification from the last completed calculation
        """
        self.report('workchain completed after {} iterations'.format(self.ctx.iteration))

        for name, port in self.spec().outputs.iteritems():
            if port.required and name not in self.ctx.restart_calc.out:
                self.report('the spec specifies the output {} as required but was not an output of {}<{}>'
                    .format(name, self.ctx.calc_name, self.ctx.restart_calc.pk))

            if name in self.ctx.restart_calc.out:
                node = self.ctx.restart_calc.out[name]
                self.out(name, self.ctx.restart_calc.out[name])
                if self._verbose:
                    self.report("attaching the node {}<{}> as '{}'"
                        .format(node.__class__.__name__, node.pk, name))

        return

    def on_terminated(self):
        """
        If the clean_workdir input was set to True, recursively collect all called Calculations by
        ourselves and our called descendants, and clean the remote folder for the JobCalculation instances
        """
        super(BaseRestartWorkChain, self).on_terminated()

        if self.inputs.clean_workdir.value is False:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = []

        for called_descendant in self.calc.called_descendants:
            if isinstance(called_descendant, JobCalculation):
                try:
                    called_descendant.out.remote_folder._clean()
                    cleaned_calcs.append(called_descendant.pk)
                except (IOError, OSError, KeyError):
                    pass

        if cleaned_calcs:
            self.report('cleaned remote folders of calculations: {}'.format(' '.join(map(str, cleaned_calcs))))

    def _handle_calculation_sanity_checks(self, calculation):
        """
        Calculations that run successfully may still have problems that can only be determined when inspecting
        the output. The same problems may also be the hidden root of a calculation failure. For that reason,
        after verifying that the calculation ran, regardless of its calculation state, we perform some sanity
        checks.
        """
        pass

    def _handle_submission_failure(self, calculation):
        """
        The submission of the calculation has failed. If the submission_failure flag is set to true, this
        is the second consecutive submission failure and we abort the workchain Otherwise we restart once more.
        """
        if self.ctx.submission_failure:
            self.report('submission for {}<{}> failed for the second consecutive time'
                .format(self.ctx.calc_name, calculation.pk))
            return self.exit_codes.ERROR_SECOND_CONSECUTIVE_SUBMISSION_FAILURE

        else:
            self.report('submission for {}<{}> failed, restarting once more'
                .format(self.ctx.calc_name, calculation.pk))

    def _handle_unexpected_failure(self, calculation, exception=None):
        """
        The calculation has failed for an unknown reason and could not be handled. If the unexpected_failure
        flag is true, this is the second consecutive unexpected failure and we abort the workchain.
        Otherwise we restart once more.
        """
        if exception:
            self.report('{}'.format(exception))

        if self.ctx.unexpected_failure:
            self.report('failure of {}<{}> could not be handled for the second consecutive time'
                .format(self.ctx.calc_name, calculation.pk))
            return self.exit_codes.ERROR_SECOND_CONSECUTIVE_UNHANDLED_FAILURE

        else:
            self.report('failure of {}<{}> could not be handled, restarting once more'
                .format(self.ctx.calc_name, calculation.pk))

    def _handle_calculation_failure(self, calculation):
        """
        The calculation has failed so we try to analyze the reason and change the inputs accordingly
        for the next calculation. If the calculation failed, but did so cleanly, we set it as the
        restart_calc, in all other cases we do not replace the restart_calc
        """
        try:
            outputs = calculation.out.output_parameters.get_dict()
            _ = outputs['warnings']
            _ = outputs['parser_warnings']
        except (AttributeError, KeyError) as exception:
            raise UnexpectedCalculationFailure(exception)

        is_handled = False

        # Sort the handlers based on their priority in reverse order
        handlers = sorted(self._error_handlers, key=lambda x: x.priority, reverse=True)

        if not handlers:
            raise UnexpectedCalculationFailure('no calculation error handlers were registered')

        for handler in handlers:
            handler_report = handler.method(self, calculation)

            # If at least one error is handled, we consider the calculation failure handled
            if handler_report and handler_report.is_handled:
                self.ctx.restart_calc = calculation
                is_handled = True

            # After certain error handlers, we may want to skip all other error handling
            if handler_report and handler_report.do_break:
                break

        # If none of the executed error handlers reported that they handled an error, the failure reason is unknown
        if not is_handled:
            raise UnexpectedCalculationFailure('calculation failure was not handled')

        return handler_report.exit_code

    def _prepare_process_inputs(self, process, inputs):
        """
        Prepare the inputs dictionary for a calculation process. Any remaining bare dictionaries in the inputs
        dictionary will be wrapped in a ParameterData data node except for the 'options' key which should remain
        a standard dictionary. Another exception are dictionaries whose keys are not strings but for example tuples.
        This is the format used by input groups as in for example the explicit pseudo dictionary where the key is
        a tuple of kind to which the UpfData corresponds.
        """
        from aiida_quantumespresso.utils.mapping import prepare_process_inputs
        return prepare_process_inputs(process, inputs)
