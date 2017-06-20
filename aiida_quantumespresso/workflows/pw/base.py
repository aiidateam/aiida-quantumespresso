# -*- coding: utf-8 -*-
from aiida.orm import Code
from aiida.orm.data.base import Bool, Int, Str
from aiida.orm.data.upf import UpfData, get_pseudos_from_structure
from aiida.orm.data.remote import RemoteData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.common.datastructures import calc_states
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, while_
from aiida_quantumespresso.calculations.pw import PwCalculation

class PwBaseWorkChain(WorkChain):
    """
    Base Workchain to launch a Quantum Espresso pw.x total energy calculation
    """
    def __init__(self, *args, **kwargs):
        super(PwBaseWorkChain, self).__init__(*args, **kwargs)

    @classmethod
    def define(cls, spec):
        super(PwBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        spec.input_group('pseudos', required=False)
        spec.input('pseudo_family', valid_type=Str, required=False)
        spec.input('parent_folder', valid_type=RemoteData, required=False)
        spec.input('kpoints', valid_type=KpointsData)
        spec.input('parameters', valid_type=ParameterData)
        spec.input('settings', valid_type=ParameterData)
        spec.input('options', valid_type=ParameterData)
        spec.input('clean_workdir', valid_type=Bool, default=Bool(False))
        spec.input('max_iterations', valid_type=Int, default=Int(10))
        spec.outline(
            cls.setup,
            cls.validate_pseudo_potentials,
            while_(cls.should_run_pw)(
                cls.run_pw,
                cls.inspect_pw,
            ),
            cls.run_results,
        )
        spec.dynamic_output()

    def setup(self):
        """
        Initialize context variables
        """
        self.ctx.max_iterations = self.inputs.max_iterations.value
        self.ctx.restart_calc = None
        self.ctx.is_finished = False
        self.ctx.iteration = 0

        # Define convenience dictionary of inputs for PwCalculation
        self.ctx.inputs = {
            'code': self.inputs.code,
            'structure': self.inputs.structure,
            'pseudo': {},
            'kpoints': self.inputs.kpoints,
            'parameters': self.inputs.parameters.get_dict(),
            'settings': self.inputs.settings.get_dict(),
            '_options': self.inputs.options.get_dict(),
        }

        # Prevent PwCalculation from being terminated by scheduler
        max_wallclock_seconds = self.ctx.inputs['_options']['max_wallclock_seconds']
        self.ctx.inputs['parameters']['CONTROL']['max_seconds'] = int(0.95 * max_wallclock_seconds)

        return

    def validate_pseudo_potentials(self):
        """
        Validate the inputs related to pseudopotentials to check that we have the minimum required
        amount of information to be able to run a PwCalculation
        """
        structure = self.inputs.structure
        pseudo_family = self.inputs.pseudo_family.value

        if all([key not in self.inputs for key in ['pseudos', 'pseudo_family']]):
            self.abort_nowait('neither explicit pseudos nor a pseudo_family was specified in the inputs')
            return
        elif all([key in self.inputs for key in ['pseudos', 'pseudo_family']]):
            self.report('both explicit pseudos as well as a pseudo_family were specified: using explicit pseudos')
            self.ctx.inputs['pseudo'] = self.inputs.pseudos
        elif 'pseudos' in self.inputs:
            self.report('only explicit pseudos were specified: using explicit pseudos')
            self.ctx.inputs['pseudo'] = self.inputs.pseudos
        elif 'pseudo_family' in self.inputs:
            self.report('only a pseudo_family was specified: using pseudos from pseudo_family')
            self.ctx.inputs['pseudo'] = get_pseudos_from_structure(structure, pseudo_family)

        for kind in self.inputs.structure.get_kind_names():
            if kind not in self.ctx.inputs['pseudo']:
                self.abort_nowait('no pseudo available for element {}'.format(kind))
            elif not isinstance(self.ctx.inputs['pseudo'][kind], UpfData):
                self.abort_nowait('pseudo for element {} is not of type UpfData'.format(kind))

    def should_run_pw(self):
        """
        Return whether a pw restart calculation should be run, which is the case as long as the last
        calculation was not converged successfully and the maximum number of restarts has not yet
        been exceeded
        """
        return not self.ctx.is_finished and self.ctx.iteration < self.ctx.max_iterations

    def run_pw(self):
        """
        Run a new PwCalculation or restart from a previous PwCalculation run in this workchain
        """
        self.ctx.iteration += 1

        # Create local copy of general inputs stored in the context and adapt for next calculation
        inputs = dict(self.ctx.inputs)

        if self.ctx.iteration == 1 and 'parent_folder' in self.inputs:
            inputs['parameters']['CONTROL']['restart_mode'] = 'restart'
            inputs['parent_folder'] = self.inputs.parent_folder
        elif self.ctx.restart_calc:
            inputs['parameters']['CONTROL']['restart_mode'] = 'restart'
            inputs['parent_folder'] = self.ctx.restart_calc.out.remote_folder
        else:
            inputs['parameters']['CONTROL']['restart_mode'] = 'from_scratch'

        inputs['parameters'] = ParameterData(dict=inputs['parameters'])
        inputs['settings'] = ParameterData(dict=inputs['settings'])

        process = PwCalculation.process()
        running = submit(process, **inputs)

        self.report('launching PwCalculation<{}> iteration #{}'.format(running.pid, self.ctx.iteration))

        return ToContext(calculation=running)

    def inspect_pw(self):
        """
        Analyse the results of the previous PwCalculation, checking whether it finished successfully
        or if not troubleshoot the cause and adapt the input parameters accordingly before
        restarting, or abort if unrecoverable error was found
        """
        try:
            calculation = self.ctx.calculation
        except Exception:
            self.abort_nowait('the first iteration finished without returning a PwCalculation')
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
            self.abort_nowait('last ran PwCalculation<{}>'.format(calculation.pk))

        # Abort: unexpected state of last calculation
        elif calculation.get_state() not in expected_states:
            self.abort_nowait('unexpected state ({}) of PwCalculation<{}>'.format(
                calculation.get_state(), calculation.pk))

        # Retry: submission failed, try to restart or abort
        elif calculation.get_state() in [calc_states.SUBMISSIONFAILED]:
            self._handle_submission_failure(calculation)

        # Retry: calculation failed, try to salvage or abort
        elif calculation.get_state() in [calc_states.FAILED]:
            self._handle_calculation_failure(calculation)

        # Retry: try to convergence restarting from this calculation
        else:
            self.report('calculation did not converge after {} iterations, restarting'.format(self.ctx.iteration))
            self.ctx.restart_calc = calculation

        return

    def run_results(self):
        """
        Attach the output parameters and retrieved folder of the last calculation to the outputs
        """
        self.report('workchain completed after {} iterations'.format(self.ctx.iteration))
        self.out('output_parameters', self.ctx.restart_calc.out.output_parameters)
        self.out('remote_folder', self.ctx.restart_calc.out.remote_folder)
        self.out('retrieved', self.ctx.restart_calc.out.retrieved)

        if 'output_structure' in self.ctx.restart_calc.out:
            self.out('output_structure', self.ctx.restart_calc.out.output_structure)

    def _handle_submission_failure(self, calculation):
        """
        The submission of the calculation has failed, if it was the second consecutive failure we
        abort the workchain, else we set the has_submission_failed flag and try again
        """
        self.abort_nowait('submission failed for the {} in iteration {}, but error handling is not implemented yet'
            .format(PwCalculation.__name__, self.ctx.iteration))

    def _handle_calculation_failure(self, calculation):
        """
        The calculation has failed so we try to analyze the reason and change the inputs accordingly
        for the next calculation. If the calculation failed, but did so cleanly, we set it as the
        restart_calc, in all other cases we do not replace the restart_calc
        """
        self.abort_nowait('execution failed for the {} in iteration {}, but error handling is not implemented yet'
            .format(PwCalculation.__name__, self.ctx.iteration))

    def on_stop(self):
        """
        Clean remote folders of the PwCalculations that were run if the clean_workdir parameter was
        set to true in the Workchain inputs
        """
        super(PwBaseWorkChain, self).on_stop()

        if not self.inputs.clean_workdir.value:
            self.report('remote folders will not be cleaned')
            return

        for calc in self.ctx.calculation:
            try:
                calc.out.remote_folder._clean()
                self.report('cleaned remote folder of {}<{}>'.format(calc.__class__.__name__, calc.pk))
            except Exception:
                pass