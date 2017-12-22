# -*- coding: utf-8 -*-
from collections import namedtuple
from aiida.common.datastructures import calc_states
from aiida.orm.calculation import JobCalculation
from aiida.orm.data.parameter import ParameterData
from aiida.work.workchain import WorkChain, ToContext, append_
from aiida.work.run import submit
from aiida_quantumespresso.common.pluginloader import get_plugin, get_plugins
from aiida_quantumespresso.common.exceptions import UnexpectedFailure

class BaseWorkChain(WorkChain):
    """
    Base Workchain

    The idea is to subclass this workchain and leverage the generic error handling that is implemented
    in the few outline methods. Part of the suggested outline would look something like the following:

        while_(cls.should_run_calculation)(
            cls.run_calculation,
            cls.inspect_calculation,
        )

    Each of these methods can of course be overriden but they should be general enough to fit most
    calculation cycles. The run_calculation method will take the inputs for the calculation process
    from the context under the key 'inputs'. The user should therefore make sure that before the
    run_calculation method is called, that the to be used inputs are stored under self.ctx.inputs.
    One can update the inputs based on the results from a prior calculation by calling an outline
    method just before the run_calculation step, for example:

        while_(cls.should_run_calculation)(
            cls.prepare_calculation,
            cls.run_calculation,
            cls.inspect_calculation,
        )

    Where in the prepare_calculation method, the inputs dictionary at self.ctx.inputs is updated
    before the next calculation will be run with those inputs.
    """
    _calculation_class = None
    _error_handler_entry_point = None

    ErrorHandlingReport = namedtuple('ErrorHandlingReport', 'is_handled do_break')

    def __init__(self, *args, **kwargs):
        super(BaseWorkChain, self).__init__(*args, **kwargs)

        if self._calculation_class is None or not issubclass(self._calculation_class, JobCalculation):
            raise ValueError('no valid JobCalculation class defined for _calculation_class attribute')

        return

    def should_run_calculation(self):
        """
        Return whether a new calculation should be run, which is the case as long as the last
        calculation was not finished successfully and the maximum number of restarts has not yet
        been exceeded
        """
        return not self.ctx.is_finished and self.ctx.iteration < self.ctx.max_iterations

    def run_calculation(self):
        """
        Run a new calculation
        """
        self.ctx.iteration += 1

        inputs = self._prepare_process_inputs(self.ctx.inputs)
        process = self._calculation_class.process()
        running = submit(process, **inputs)

        self.report('launching {}<{}> iteration #{}'
            .format(self._calculation_class.__name__, running.pid, self.ctx.iteration))

        return ToContext(calculations=append_(running))

    def inspect_calculation(self):
        """
        Analyse the results of the previous calculation, checking whether it finished successfully
        or if not troubleshoot the cause and handle the errors, or abort if unrecoverable error was found
        """
        try:
            calculation = self.ctx.calculations[-1]
        except IndexError:
            self.abort_nowait('the first iteration finished without returning a {}'
                .format(self._calculation_class.__name__))
            return

        expected_states = [calc_states.FINISHED, calc_states.FAILED, calc_states.SUBMISSIONFAILED]

        # Done: successful convergence of last calculation
        if calculation.has_finished_ok():
            self.report('converged successfully after {} iterations'.format(self.ctx.iteration))
            self.ctx.restart_calc = calculation
            self.ctx.is_finished = True

        # Abort: exceeded maximum number of retries
        elif self.ctx.iteration >= self.ctx.max_iterations:
            self.report('reached the maximumm number of iterations {}'.format(self.ctx.max_iterations))
            self.abort_nowait('last ran {}<{}>'.format(self._calculation_class.__name__, calculation.pk))

        # Abort: unexpected state of last calculation
        elif calculation.get_state() not in expected_states:
            self.abort_nowait('unexpected state ({}) of {}<{}>'
                .format(calculation.get_state(), self._calculation_class.__name__, calculation.pk))

        # Retry or abort: submission failed, try to restart or abort
        elif calculation.get_state() in [calc_states.SUBMISSIONFAILED]:
            self._handle_submission_failure(calculation)
            self.ctx.submission_failure = True

        # Retry or abort: calculation finished or failed
        elif calculation.get_state() in [calc_states.FINISHED, calc_states.FAILED]:

            # Calculation was at least submitted successfully, so we reset the flag
            self.ctx.submission_failure = False

            # Check output for problems independent on calculation state and that do not trigger parser warnings
            self._handle_calculation_sanity_checks(calculation)

            # Calculation failed, try to salvage it or handle any unexpected failures
            if calculation.get_state() in [calc_states.FAILED]:
                try:
                    handled = self._handle_calculation_failure(calculation)
                except UnexpectedFailure as exception:
                    self._handle_unexpected_failure(calculation, exception)
                    self.ctx.unexpected_failure = True

            # Calculation finished: but did not finish ok, simply try to restart from this calculation
            else:
                self.ctx.unexpected_failure = False
                self.ctx.restart_calc = calculation
                self.report('calculation did not converge after {} iterations, restarting'
                    .format(self.ctx.iteration))

        return

    def results(self):
        """
        Attach the outputs specified in the output specification from the last completed calculation
        """
        self.report('workchain completed after {} iterations'.format(self.ctx.iteration))

        for name, port in self.spec().outputs.iteritems():
            if port.required and name not in self.ctx.restart_calc.out:
                self.report('the spec specifies the output {} as required but was not an output of {}<{}>'
                    .format(name, self._calculation_class.__name__, self.ctx.restart_calc.pk))

            if name in self.ctx.restart_calc.out:
                node = self.ctx.restart_calc.out[name]
                self.out(name, self.ctx.restart_calc.out[name])
                self.report("attaching the node {}<{}> as '{}'"
                    .format(node.__class__.__name__, node.pk, name))

        return

    def on_destroy(self):
        """
        Clean remote folders of the calculations that were run if the clean_workdir input was
        provided and True in the Workchain inputs
        """
        super(BaseWorkChain, self).on_destroy()
        if not self.has_finished():
            return

        if not 'clean_workdir' in self.inputs or self.inputs.clean_workdir.value is False:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = []

        try:
            calculations = self.ctx.calculations
        except AttributeError:
            calculations = []

        for calculation in calculations:
            try:
                calculation.out.remote_folder._clean()
                cleaned_calcs.append(calculation.pk)
            except Exception:
                pass

        if len(cleaned_calcs) > 0:
            self.report('cleaned remote folders of calculations: {}'.format(' '.join(map(str, cleaned_calcs))))

    def _handle_calculation_sanity_checks(self, calculation):
        """
        Calculations that run successfully may still have problems that can only be determined when inspecting
        the output. The same problems may also be the hidden root of a calculation failure. For that reason,
        after verifying that the calculation ran, regardless of its calculation state, we perform some sanity
        checks.
        """
        return

    def _handle_submission_failure(self, calculation):
        """
        The submission of the calculation has failed. If the submission_failure flag is set to true, this
        is the second consecutive submission failure and we abort the workchain Otherwise we restart once more.
        """
        if self.ctx.submission_failure:
            self.abort_nowait('submission for {}<{}> failed for the second consecutive time'
                .format(calculation.pk))
        else:
            self.report('submission for {}<{}> failed, restarting once more'
                .format(self._calculation_class.__name__, calculation.pk))

        return

    def _handle_unexpected_failure(self, calculation, exception=None):
        """
        The calculation has failed for an unknown reason and could not be handled. If the unexpected_failure
        flag is true, this is the second consecutive unexpected failure and we abort the workchain.
        Otherwise we restart once more.
        """
        if exception:
            self.report('exception: {}'.format(exception))

        if self.ctx.unexpected_failure:
            self.abort_nowait('{}<{}> failed for an unknown case for the second consecutive time'
                .format(self._calculation_class.__name__, calculation.pk))
        else:
            self.report('{}<{}> failed for an unknown case, restarting once more'
                .format(self._calculation_class.__name__, calculation.pk))

        return

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
        except (NotExistent, AttributeError, KeyError) as exception:
            raise UnexpectedFailure(exception)

        is_handled = False

        if self._error_handler_entry_point is None:
            self.report('no error handler entry point registered cannot handle calculation failure')
            return

        error_handlers = []
        for plugin in get_plugins(self._error_handler_entry_point):
            plugin_error_handlers = get_plugin(self._error_handler_entry_point, plugin)()
            error_handlers.extend(plugin_error_handlers)

        for handler in error_handlers:
            handler_report = handler(self, calculation)

            # If at least one error is handled, we consider the computation failure handled
            if handler_report and handler_report.is_handled:
                self.ctx.restart_calc = calculation
                is_handled = True

            # After certain error handlers, we may want to skip all other error handling
            if handler_report and handler_report.do_break:
                break

        # If none of the executed error handler reported that they handled an error, the failure reason is unknown
        if not is_handled:
            raise UnexpectedFailure('{}<{}> failed due to an unknown reason'
                .format(self._calculation_class.__name__, calculation.pk))

        return

    def _prepare_process_inputs(self, inputs=None):
        """
        Prepare the inputs dictionary for a calculation process. The 'max_seconds' setting in the 'CONTROL' card
        of the parameters will be set to a fraction of the 'max_wallclock_seconds' that will be given to the job via
        the '_options' dictionary. This will prevent the job from being prematurely terminated by the scheduler without
        getting the chance to exit cleanly. Any remaining bare dictionaries in the inputs dictionary will be wrapped
        in a ParameterData data node except for the '_options' key which should remain a standard dictionary
        """
        if inputs == None:
            inputs = self.ctx.inputs

        prepared_inputs = {}

        # Limit the max seconds to a fraction of the scheduler's max_wallclock_seconds to prevent early termination
        max_wallclock_seconds = inputs['_options']['max_wallclock_seconds']
        max_seconds_factor = self.defaults.delta_factor_max_seconds
        max_seconds = max_wallclock_seconds * max_seconds_factor
        inputs['parameters']['CONTROL']['max_seconds'] = max_seconds

        # Wrap all the bare dictionaries in a ParameterData
        for key, value in inputs.iteritems():
            if key != '_options' and isinstance(value, dict) and all([isinstance(k, (str, unicode)) for k in value.keys()]):
                prepared_inputs[key] = ParameterData(dict=value)
            else:
                prepared_inputs[key] = value

        return prepared_inputs