# -*- coding: utf-8 -*-
from aiida.orm import Code
from aiida.orm.data.base import Bool, Float, Int, Str
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.utils import CalculationFactory
from aiida.common.exceptions import NotExistent
from aiida.common.datastructures import calc_states
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, while_, append_
from aiida_quantumespresso.common.exceptions import UnexpectedFailure

PhCalculation = CalculationFactory('quantumespresso.ph')
PwCalculation = CalculationFactory('quantumespresso.pw')

class PhBaseWorkChain(WorkChain):
    """
    Base Workchain to launch a Quantum Espresso phonon ph.x calculation and restart it until
    successfully converged or until the maximum number of restarts is exceeded
    """
    def __init__(self, *args, **kwargs):
        super(PhBaseWorkChain, self).__init__(*args, **kwargs)

    @classmethod
    def define(cls, spec):
        super(PhBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('parent_calc', valid_type=PwCalculation)
        spec.input('qpoints', valid_type=KpointsData)
        spec.input('parameters', valid_type=ParameterData)
        spec.input('settings', valid_type=ParameterData)
        spec.input('options', valid_type=ParameterData)
        spec.input('clean_workdir', valid_type=Bool, default=Bool(False))
        spec.input('compute_epsil', valid_type=Bool, default=Bool(False))
        spec.input('alpha_mix', valid_type=Float, default=Float(0.7))
        spec.input('max_iterations', valid_type=Int, default=Int(10))
        spec.outline(
            cls.setup,
            while_(cls.should_run_ph)(
                cls.run_ph,
                cls.inspect_ph,
            ),
            cls.run_results,
            cls.run_clean,
        )
        spec.dynamic_output()

    def setup(self):
        """
        Initialize context variables
        """
        self.ctx.max_iterations = self.inputs.max_iterations.value
        self.ctx.restart_calc = self.inputs.parent_calc
        self.ctx.unexpected_failure = False
        self.ctx.submission_failure = False
        self.ctx.is_finished = False
        self.ctx.iteration = 0

        # Define convenience dictionary of inputs for PhCalculation
        self.ctx.inputs = {
            'code': self.inputs.code,
            'qpoints': self.inputs.qpoints,
            'parameters': self.inputs.parameters.get_dict(),
            'settings': self.inputs.settings.get_dict(),
            '_options': self.inputs.options.get_dict(),
        }

        # Prevent PhCalculation from being terminated by scheduler
        max_wallclock_seconds = self.ctx.inputs['_options']['max_wallclock_seconds']
        self.ctx.inputs['parameters']['INPUTPH']['max_seconds'] = int(0.95 * max_wallclock_seconds)

        return

    def should_run_ph(self):
        """
        Return whether a ph restart calculation should be run, which is the case as long as the last
        calculation was not converged successfully and the maximum number of restarts has not yet
        been exceeded
        """
        return not self.ctx.is_finished and self.ctx.iteration < self.ctx.max_iterations

    def run_ph(self):
        """
        Run a PhCalculation either starting from a previous PwCalculation or restarting from a
        previous PhCalculation run in this workchain
        """
        self.ctx.iteration += 1

        # Create local copy of general inputs stored in the context and adapt for next calculation
        inputs = dict(self.ctx.inputs)

        if isinstance(self.ctx.restart_calc, PhCalculation):
            inputs['parameters']['INPUTPH']['recover'] = True
        elif isinstance(self.ctx.restart_calc, PwCalculation):
            pass
        else:
            ctype = type(self.ctx.restart_calc)
            self.abort_nowait("The type '{}' of the parent calculation is invalid".format(ctype))

        inputs['parent_folder'] = self.ctx.restart_calc.out.remote_folder
        inputs['parameters'] = ParameterData(dict=inputs['parameters'])
        inputs['settings'] = ParameterData(dict=inputs['settings'])

        process = PhCalculation.process()
        running = submit(process, **inputs)

        self.report('launching PhCalculation<{}> iteration #{}'.format(
            running.pid, self.ctx.iteration))

        return ToContext(calculations=append_(running))

    def inspect_ph(self):
        """
        Analyse the results of the previous PhCalculation, checking whether it finished successfully
        or if not troubleshoot the cause and adapt the input parameters accordingly before
        restarting, or abort if unrecoverable error was found
        """
        try:
            calculation = self.ctx.calculations[-1]
        except IndexError:
            self.abort_nowait('the first iteration finished without returning a PhCalculation')
            return

        expected_states = [calc_states.FINISHED, calc_states.FAILED, calc_states.SUBMISSIONFAILED]

        # Done: successful convergence of last calculation
        if calculation.has_finished_ok():
            self.report('converged successfully after {} iterations'.format(self.ctx.iteration))
            self.ctx.restart_calc = calculation
            self.ctx.is_finished = True

        # Abort: exceeded maximum number of retries
        elif self.ctx.iteration >= self.ctx.max_iterations:
            self.report('reached the max number of iterations {}'.format(self.ctx.max_iterations))
            self.abort_nowait('last ran PhCalculation<{}>'.format(calculation.pk))

        # Abort: unexpected state of last calculation
        elif calculation.get_state() not in expected_states:
            self.abort_nowait('unexpected state ({}) of PhCalculation<{}>'.format(
                calculation.get_state(), calculation.pk))

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
                    self._handle_calculation_failure(calculation)
                except UnexpectedFailure as exception:
                    self._handle_unexpected_failure(calculation, exception)
                    self.ctx.unexpected_failure = True

            # Calculation finished: but did not finish ok, simply try to restart from this calculation
            else:
                self.ctx.unexpected_failure = False
                self.ctx.restart_calc = calculation
                self.report('calculation did not converge after {} iterations, restarting'.format(self.ctx.iteration))

        return

    def run_results(self):
        """
        Attach the output parameters and retrieved folder of the last calculation to the outputs
        """
        calculation = self.ctx.calculations[-1]
        self.report('workchain completed after {} iterations'.format(self.ctx.iteration))
        self.out('output_parameters', calculation.out.output_parameters)
        self.out('retrieved', calculation.out.retrieved)

    def run_clean(self):
        """
        Clean remote folders of the PhCalculations that were run if the clean_workdir parameter was
        set to true in the Workchain inputs
        """
        if not self.inputs.clean_workdir.value:
            self.report('remote folders will not be cleaned')
            return

        for calc in self.ctx.calculations:
            try:
                calc.out.remote_folder._clean()
                self.report('cleaned remote folder of {}<{}>'.format(calc.__class__.__name__, calc.pk))
            except Exception:
                pass

    def _handle_submission_failure(self, calculation):
        """
        The submission of the calculation has failed. If the submission_failure flag is set to true, this
        is the second consecutive submission failure and we abort the workchain. Otherwise we restart once more.
        """
        if self.ctx.submission_failure:
            self.abort_nowait('submission for PhCalculation<{}> failed for the second consecutive time'
                .format(calculation.pk))
        else:
            self.report('submission for PhCalculation<{}> failed, restarting once more'
                .format(calculation.pk))

    def _handle_unexpected_failure(self, calculation, exception=None):
        """
        The calculation has failed for an unknown reason and could not be handled. If the unexpected_failure
        flag is true, this is the second consecutive unexpected failure and we abort the workchain.
        Otherwise we restart once more.
        """
        if exception:
            self.report('exception: {}'.format(exception))

        if self.ctx.unexpected_failure:
            self.abort_nowait('PhCalculation<{}> failed for an unknown case for the second consecutive time'
                .format(calculation.pk))
        else:
            self.report('PhCalculation<{}> failed for an unknown case, restarting once more'
                .format(calculation.pk))

    def _handle_calculation_sanity_checks(self, calculation):
        """
        Calculations that run successfully may still have problems that can only be determined when inspecting
        the output. The same problems may also be the hidden root of a calculation failure. For that reason,
        after verifying that the calculation ran, regardless of its calculation state, we perform some sanity
        checks.
        """
        return

    def _handle_calculation_failure(self, calculation):
        """
        The calculation has failed so we try to analyze the reason and change the inputs accordingly
        for the next calculation. If the calculation failed, but did so cleanly, we set it as the
        restart_calc, in all other cases we do not replace the restart_calc
        """
        try:
            input_parameters = calculation.inp.parameters.get_dict()
            output_parameters = calculation.out.output_parameters.get_dict()
            warnings = output_parameters['warnings']
            parser_warnings = output_parameters['parser_warnings']
        except (NotExistent, AttributeError, KeyError) as exception:
            raise UnexpectedFailure(exception)
        else:

            # Errors during parsing of input files
            if any(['read_namelists' in w for w in warnings]):
                self._handle_fatal_error_read_namelists(calculation)

            # Premature termination by the scheduler
            elif (('QE ph run did not reach the end of the execution.' in parser_warnings) and len(warnings) == 0):
                self._handle_error_premature_termination(calculation)

            # Calculation did not converge
            elif ('Phonon did not reach end of self consistency' in warnings):
                self._handle_fatal_error_not_converged(calculation)

            # Clean calculation termination due to exceeding maximum allotted walltime
            elif 'Maximum CPU time exceeded' in warnings:
                self._handle_error_exceeded_maximum_walltime(calculation)

            else:
                raise UnexpectedFailure

    def _handle_fatal_error_read_namelists(self, calculation):
        """
        The calculation failed because it could not read the generated input file
        """
        self.abort_nowait('PhCalculation<{}> failed because of an invalid input file'
            .format(calculation.pk))

    def _handle_fatal_error_not_converged(self, calculation):
        """
        The calculation failed because it could not read the generated input file
        """
        alpha_mix_old = calculation.inp.parameters.get_dict()['INPUTPH'].get('alpha_mix(1)', self.inputs.alpha_mix.value)
        alpha_mix_new = 0.9 * alpha_mix_old
        self.ctx.inputs['parameters']['INPUTPH']['alpha_mix(1)'] = alpha_mix_new
        self.ctx.restart_calc = calculation
        self.report('PhCalculation<{}> terminated without reaching convergence, '
            'setting alpha_mix to {} and restarting'.format(calculation.pk, alpha_mix_new))

    def _handle_error_premature_termination(self, calculation):
        """
        Calculation did not reach the end of execution, probably because it was killed by the scheduler
        for running out of allotted walltime
        """
        input_parameters = calculation.inp.parameters.get_dict()
        settings = self.ctx.inputs['settings']
        max_seconds = settings.get('max_seconds', input_parameters['CONTROL']['max_seconds'])
        max_seconds_reduced = int(0.95 * max_seconds)
        self.ctx.inputs['parameters']['CONTROL']['max_seconds'] = max_seconds_reduced
        self.report('PhCalculation<{}> was terminated prematurely, reducing "max_seconds" from {} to {}'
            .format(calculation.pk, max_seconds, max_seconds_reduced))

    def _handle_error_exceeded_maximum_walltime(self, calculation):
        """
        Calculation ended nominally but ran out of allotted wall time
        """
        self.ctx.restart_calc = calculation
        self.report('PhCalculation<{}> terminated because maximum wall time was exceeded, restarting'
            .format(calculation.pk))