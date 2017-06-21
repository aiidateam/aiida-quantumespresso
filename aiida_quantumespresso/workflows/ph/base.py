# -*- coding: utf-8 -*-
from aiida.orm import Code
from aiida.orm.data.base import Bool, Float, Int, Str
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.common.datastructures import calc_states
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, while_, append_
from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.calculations.pw import PwCalculation

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
        )
        spec.dynamic_output()

    def setup(self):
        """
        Initialize context variables
        """
        self.ctx.max_iterations = self.inputs.max_iterations.value
        self.ctx.restart_calc = self.inputs.parent_calc
        self.ctx.has_calculation_failed = False
        self.ctx.has_submission_failed = False
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
        except Exception:
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

        # Retry: submission failed, try to restart or abort
        elif calculation.get_state() in [calc_states.SUBMISSIONFAILED]:
            self._handle_submission_failure(calculation)
            self.ctx.has_submission_failed = True

        # Retry: calculation failed, try to salvage or abort
        elif calculation.get_state() in [calc_states.FAILED]:
            self._handle_calculation_failure(calculation)
            self.ctx.has_submission_failed = False

        return

    def run_results(self):
        """
        Attach the output parameters and retrieved folder of the last calculation to the outputs
        """
        calculation = self.ctx.calculations[-1]
        self.report('workchain completed after {} iterations'.format(self.ctx.iteration))
        self.out('output_parameters', calculation.out.output_parameters)
        self.out('retrieved', calculation.out.retrieved)

    def on_stop(self):
        """
        Clean remote folders of the PhCalculations that were run if the clean_workdir parameter was
        set to true in the Workchain inputs
        """
        super(PhBaseWorkChain, self).on_stop()

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
        The submission of the calculation has failed, if it was the second consecutive failure we
        abort the workchain, else we set the has_submission_failed flag and try again
        """
        if self.ctx.has_submission_failed:
            self.abort_nowait('submission for PhCalculation<{}> failed for the second time'.format(
                calculation.pk))
        else:
            self.report('submission for PhCalculation<{}> failed, retrying once more'.format(
                calculation.pk))

    def _handle_calculation_failure(self, calculation):
        """
        The calculation has failed so we try to analyze the reason and change the inputs accordingly
        for the next calculation. If the calculation failed, but did so cleanly, we set it as the
        restart_calc, in all other cases we do not replace the restart_calc
        """
        if any(['namelist' in w for w in calculation.res.warnings]):
            self.ctx.has_calculation_failed = False
            warnings = '||'.join([w.strip() for w in calculation.res.warnings])
            self.report('PhCalculation<{}> failed due to incorrect input file'.format(calculation.pk))
            self.report('list of warnings: {}'.format(warnings))
            self.abort_nowait('workchain failed after {} iterations'.format(self.ctx.iteration))

        elif (('QE ph run did not reach the end of the execution.' in calculation.res.parser_warnings)
            and len(calculation.res.warnings) == 0):
            self.ctx.has_calculation_failed = False
            max_seconds_old = calculation.inp.parameters.get_dict()['INPUTPH']['max_seconds']
            max_seconds_new = int(0.95 * max_seconds_old)
            self.ctx.inputs['parameters']['INPUTPH']['max_seconds'] = max_seconds_new
            self.report('PhCalculation<{}> terminated in the middle of scf step, '
                'setting max_seconds to {}'.format(calculation.pk, max_seconds_new))

        elif ('Phonon did not reach end of self consistency' in calculation.res.warnings):
            self.ctx.has_calculation_failed = False
            alpha_mix_old = calculation.inp.parameters.get_dict()['INPUTPH'].get('alpha_mix(1)',
                self.inputs.alpha_mix.value)
            alpha_mix_new = 0.9 * alpha_mix_old
            self.ctx.inputs['parameters']['INPUTPH']['alpha_mix(1)'] = alpha_mix_new
            self.ctx.restart_calc = calculation
            self.report('PhCalculation<{}> terminated without reaching convergence, '
                'setting alpha_mix to {} and restarting'.format(calculation.pk, alpha_mix_new))

        elif 'Maximum CPU time exceeded' in calculation.res.warnings:
            self.ctx.has_calculation_failed = False
            self.ctx.restart_calc = calculation
            self.report('PhCalculation<{}> terminated because maximum walltime was exceeded, '
                'restarting'.format(calculation.pk))

        else:

            # If has_calculation_failed is set, second unexpected failure in a row, so we abort
            if self.ctx.has_calculation_failed:
                warnings = '||'.join([w.strip() for w in calculation.res.warnings])
                parser_warnings = '||'.join([w.strip() for w in calculation.res.parser_warnings])
                self.report('PhCalculation<{}> failed unexpectedly'.format(calculation.pk))
                self.report('list of warnings: {}'.format(warnings))
                self.report('list of parser warnings: {}'.format(parser_warnings))
                self.abort_nowait('second unexpected failure in a row'.format(calculation.pk))
            else:
                self.ctx.has_calculation_failed = True
                self.report('PhCalculation<{}> failed unexpectedly, '
                    'restarting only once more'.format(calculation.pk))