# -*- coding: utf-8 -*-
from collections import namedtuple
from aiida.orm import Code
from aiida.orm.data.base import Bool, Int, Str
from aiida.orm.data.upf import UpfData, get_pseudos_from_structure
from aiida.orm.data.folder import FolderData
from aiida.orm.data.remote import RemoteData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.singlefile import SinglefileData
from aiida.orm.utils import CalculationFactory
from aiida.common.exceptions import NotExistent
from aiida.common.datastructures import calc_states
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, while_, append_
from aiida_quantumespresso.common.exceptions import UnexpectedFailure
from aiida_quantumespresso.common.pluginloader import get_plugin, get_plugins

PwCalculation = CalculationFactory('quantumespresso.pw')

class PwBaseWorkChain(WorkChain):
    """
    Base Workchain to launch a Quantum Espresso pw.x total energy calculation
    """

    ErrorHandlingReport = namedtuple('ErrorHandlingReport', 'is_handled do_break')

    def __init__(self, *args, **kwargs):
        super(PwBaseWorkChain, self).__init__(*args, **kwargs)

        # Default values
        self.defaults = {
            'qe': {
                'degauss': 0.,
                'diagonalization': 'david',
                'electron_maxstep': 100,
                'mixing_beta': 0.7,
                'mixing_mode': 'plain',
                'mixing_ndim': 8,
                'noncolin': False,
                'nspin': 1,
                'occupations': None,
                'press': 0.,
                'press_conv_thr': 0.5,
                'smearing': '',
                'startmag': 0.,
                'wf_collect': False,
            },
            'delta_threshold_degauss': 30,
            'delta_factor_degauss': 0.1,
            'delta_factor_mixing_beta': 0.8,
            'delta_factor_max_seconds': 0.95,
        }

    @classmethod
    def define(cls, spec):
        super(PwBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        spec.input_group('pseudos', required=False)
        spec.input('pseudo_family', valid_type=Str, required=False)
        spec.input('parent_folder', valid_type=RemoteData, required=False)
        spec.input('kpoints', valid_type=KpointsData)
        spec.input('vdw_table', valid_type=SinglefileData, required=False)
        spec.input('parameters', valid_type=ParameterData)
        spec.input('settings', valid_type=ParameterData)
        spec.input('options', valid_type=ParameterData)
        spec.input('clean_workdir', valid_type=Bool, default=Bool(False))
        spec.input('max_iterations', valid_type=Int, default=Int(5))
        spec.outline(
            cls.setup,
            cls.validate_pseudo_potentials,
            while_(cls.should_run_pw)(
                cls.run_pw,
                cls.inspect_pw,
            ),
            cls.results,
            cls.clean,
        )
        spec.output('output_structure', valid_type=StructureData, required=False)
        spec.output('output_parameters', valid_type=ParameterData)
        spec.output('remote_folder', valid_type=RemoteData)
        spec.output('retrieved', valid_type=FolderData)

    def setup(self):
        """
        Initialize context variables
        """
        self.ctx.max_iterations = self.inputs.max_iterations.value
        self.ctx.unexpected_failure = False
        self.ctx.submission_failure = False
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

        # Add the van der Waals kernel table file if specified
        if 'vdw_table' in self.inputs:
            self.ctx.inputs['vdw_table'] = self.inputs.vdw_table

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
            pseudos = self.inputs.pseudos
        elif 'pseudos' in self.inputs:
            self.report('only explicit pseudos were specified: using explicit pseudos')
            pseudos = self.inputs.pseudos
        elif 'pseudo_family' in self.inputs:
            self.report('only a pseudo_family was specified: using pseudos from pseudo_family {}'.format(pseudo_family))
            pseudos = get_pseudos_from_structure(structure, pseudo_family)

        for kind in self.inputs.structure.get_kind_names():
            if kind not in pseudos:
                self.abort_nowait('no pseudo available for element {}'.format(kind))
            elif not isinstance(pseudos[kind], UpfData):
                self.abort_nowait('pseudo for element {} is not of type UpfData'.format(kind))

        # The pseudos dictionary should now be a dictionary of UPF nodes with the kind as linkname
        # As such, if there are multiple kinds with the same element, there will be duplicate UPF nodes
        # but multiple links for the same input node are not allowed. Moreover, to couple the UPF nodes
        # to the Calculation instance, we have to go through the use_pseudo method, which takes the kind
        # name as an additional parameter. When creating a Calculation through a Process instance, one
        # cannot call the use methods directly but rather should pass them as keyword arguments. However, 
        # we can pass the additional parameters by using them as the keys of a dictionary
        unique_pseudos = {}
        for kind, pseudo in pseudos.iteritems():
            unique_pseudos.setdefault(pseudo, []).append(kind)

        for pseudo, kinds in unique_pseudos.iteritems():
             self.ctx.inputs['pseudo'][tuple(kinds)] = pseudo

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

        return ToContext(calculations=append_(running))

    def inspect_pw(self):
        """
        Analyse the results of the previous PwCalculation, checking whether it finished successfully
        or if not troubleshoot the cause and adapt the input parameters accordingly before
        restarting, or abort if unrecoverable error was found
        """
        try:
            calculation = self.ctx.calculations[-1]
        except IndexError:
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
                self.report('calculation did not converge after {} iterations, restarting'.format(self.ctx.iteration))

        return

    def results(self):
        """
        Attach the output parameters and retrieved folder of the last calculation to the outputs
        """
        self.report('workchain completed after {} iterations'.format(self.ctx.iteration))
        self.out('output_parameters', self.ctx.restart_calc.out.output_parameters)
        self.out('remote_folder', self.ctx.restart_calc.out.remote_folder)
        self.out('retrieved', self.ctx.restart_calc.out.retrieved)

        if 'output_structure' in self.ctx.restart_calc.out:
            self.out('output_structure', self.ctx.restart_calc.out.output_structure)

    def clean(self):
        """
        Clean remote folders of the PwCalculations that were run if the clean_workdir parameter was
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
            self.abort_nowait('submission for PwCalculation<{}> failed for the second consecutive time'
                .format(calculation.pk))
        else:
            self.report('submission for PwCalculation<{}> failed, restarting once more'
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
            self.abort_nowait('PwCalculation<{}> failed for an unknown case for the second consecutive time'
                .format(calculation.pk))
        else:
            self.report('PwCalculation<{}> failed for an unknown case, restarting once more'
                .format(calculation.pk))

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

        error_handlers = []
        for plugin in get_plugins('aiida_quantumespresso.workflows.error_handlers'):
            plugin_error_handlers = get_plugin('aiida_quantumespresso.workflows.error_handlers', plugin)()
            error_handlers.extend(plugin_error_handlers)

        for handler in error_handlers:
            handler_report = handler(self, calculation)

            # If at least one error is handled, we consider the computation failure handled
            if handler_report and handler_report.is_handled:
                is_handled = True

            # After certain error handlers, we may want to skip all other error handling
            if handler_report and handler_report.do_break:
                break

        # If none of the executed error handler reported that they handled an error, the failure reason is unknown
        if not is_handled:
            raise UnexpectedFailure('PwCalculation<{}> failed due to an unknown reason'.format(calculation.pk))


def get_error_handlers():
    """
    Return a list of all the implemented error handlers in the case of a PwCalculation failure
    """
    return [
        _handle_error_read_namelists,
        _handle_error_diagonalization,
        _handle_error_unrecognized_by_parser,
        _handle_error_exceeded_maximum_walltime
    ]

def _handle_error_read_namelists(workchain, calculation):
    """
    The calculation failed because it could not read the generated input file
    """
    if any(['read_namelists' in w for w in calculation.res.warnings]):
        workchain.abort_nowait('PwCalculation<{}> failed because of an invalid input file'.format(calculation.pk))
        return workchain.ErrorHandlingReport(True, False)

def _handle_error_diagonalization(workchain, calculation):
    """
    Diagonalization failed with current scheme. Try to restart from previous clean calculation with different scheme
    """
    input_parameters = calculation.inp.parameters.get_dict()
    input_electrons = input_parameters.get('ELECTRONS', {})
    diagonalization = input_electrons.get('diagonalization', workchain.defaults['qe']['diagonalization'])

    if ((
        any(['too many bands are not converged' in w for w in calculation.res.warnings]) or
        any(['eigenvalues not converged' in w for w in calculation.res.warnings])
    ) and (
        diagonalization == 'david'
    )):
        new_diagonalization = 'cg'
        workchain.ctx.inputs['parameters']['ELECTRONS']['diagonalization'] = 'cg'
        workchain.report('PwCalculation<{}> failed to diagonalize with "{}" scheme'.format(diagonalization))
        workchain.report('Restarting with diagonalization scheme "{}"'.format(new_diagonalization))
        return workchain.ErrorHandlingReport(True, False)

def _handle_error_unrecognized_by_parser(workchain, calculation):
    """
    Calculation failed with an error that was not recognized by the parser and was attached
    wholesale to the warnings. We treat it as an unexpected failure and raise the exception
    """
    warnings = calculation.res.warnings
    if (any(['%%%' in w for w in warnings]) or any(['Error' in w for w in warnings])):
        raise UnexpectedFailure('PwCalculation<{}> failed due to an unknown reason'.format(calculation.pk))

def _handle_error_exceeded_maximum_walltime(workchain, calculation):
    """
    Calculation ended nominally but ran out of allotted wall time
    """
    if 'Maximum CPU time exceeded' in calculation.res.warnings:
        workchain.ctx.restart_calc = calculation
        workchain.report('PwCalculation<{}> terminated because maximum wall time was exceeded, restarting'
            .format(calculation.pk))
        return workchain.ErrorHandlingReport(True, False)