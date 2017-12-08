# -*- coding: utf-8 -*-
from copy import deepcopy
from collections import namedtuple
from aiida.orm import Code
from aiida.orm.data.base import Bool, Int, Str
from aiida.orm.data.folder import FolderData
from aiida.orm.data.remote import RemoteData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.singlefile import SinglefileData
from aiida.orm.utils import CalculationFactory
from aiida.common.extendeddicts import AttributeDict
from aiida.common.exceptions import NotExistent
from aiida.common.datastructures import calc_states
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, if_, while_, append_
from aiida_quantumespresso.common.exceptions import UnexpectedFailure
from aiida_quantumespresso.common.pluginloader import get_plugin, get_plugins
from aiida_quantumespresso.utils.defaults.calculation import pw as qe_defaults
from aiida_quantumespresso.utils.mapping import update_mapping
from aiida_quantumespresso.utils.pseudopotential import validate_and_prepare_pseudos_inputs
from aiida_quantumespresso.utils.resources import get_default_options
from aiida_quantumespresso.utils.resources import get_pw_parallelization_parameters
from aiida_quantumespresso.utils.resources import cmdline_remove_npools
from aiida_quantumespresso.utils.resources import create_scheduler_resources

PwCalculation = CalculationFactory('quantumespresso.pw')

class PwBaseWorkChain(WorkChain):
    """
    Base Workchain to launch a Quantum Espresso pw.x total energy calculation
    """

    ErrorHandlingReport = namedtuple('ErrorHandlingReport', 'is_handled do_break')

    def __init__(self, *args, **kwargs):
        super(PwBaseWorkChain, self).__init__(*args, **kwargs)

        # Default values
        self.defaults = AttributeDict({
            'qe': qe_defaults,
            'delta_threshold_degauss': 30,
            'delta_factor_degauss': 0.1,
            'delta_factor_mixing_beta': 0.8,
            'delta_factor_max_seconds': 0.95,
        })

    @classmethod
    def define(cls, spec):
        super(PwBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        spec.input('kpoints', valid_type=KpointsData)
        spec.input('parameters', valid_type=ParameterData)
        spec.input_group('pseudos', required=False)
        spec.input('pseudo_family', valid_type=Str, required=False)
        spec.input('parent_folder', valid_type=RemoteData, required=False)
        spec.input('vdw_table', valid_type=SinglefileData, required=False)
        spec.input('settings', valid_type=ParameterData, required=False)
        spec.input('options', valid_type=ParameterData, required=False)
        spec.input('automatic_parallelization', valid_type=ParameterData, required=False)
        spec.input('max_iterations', valid_type=Int, default=Int(5))
        spec.input('clean_workdir', valid_type=Bool, default=Bool(False))
        spec.outline(
            cls.setup,
            cls.validate_inputs,
            if_(cls.should_run_init)(
                cls.validate_init_inputs,
                cls.run_init,
                cls.inspect_init,
            ),
            while_(cls.should_run_pw)(
                cls.run_pw,
                cls.inspect_pw,
            ),
            cls.results,
            cls.clean,
        )
        spec.output('output_band', valid_type=BandsData, required=False)
        spec.output('output_structure', valid_type=StructureData, required=False)
        spec.output('output_parameters', valid_type=ParameterData)
        spec.output('remote_folder', valid_type=RemoteData)
        spec.output('retrieved', valid_type=FolderData)

    def setup(self):
        """
        Initialize context variables and define convenience dictionary of inputs for PwCalculation. Only the required inputs
        are added here as the non required ones will have to be validated first in the next step of the outline. ParameterData
        nodes that may need to be update during the workchain are unpacked into their dictionary for convenience.
        """
        self.ctx.max_iterations = self.inputs.max_iterations.value
        self.ctx.unexpected_failure = False
        self.ctx.submission_failure = False
        self.ctx.restart_calc = None
        self.ctx.is_finished = False
        self.ctx.iteration = 0

        self.ctx.inputs = AttributeDict({
            'code': self.inputs.code,
            'structure': self.inputs.structure,
            'kpoints': self.inputs.kpoints,
            'parameters': self.inputs.parameters.get_dict()
        })

        return

    def validate_inputs(self):
        """
        Validate inputs that may depend on each other
        """
        if 'CONTROL'not in self.ctx.inputs.parameters:
            self.ctx.inputs.parameters['CONTROL'] = {}

        if 'calculation' not in self.ctx.inputs.parameters['CONTROL']:
            self.ctx.inputs.parameters['CONTROL']['calculation'] = 'scf'

        if 'parent_folder' in self.inputs:
            self.ctx.inputs.parent_folder = self.inputs.parent_folder
            self.ctx.inputs.parameters['CONTROL']['restart_mode'] = 'restart'
        else:
            self.ctx.inputs.parent_folder = None
            self.ctx.inputs.parameters['CONTROL']['restart_mode'] = 'from_scratch'

        if 'settings' in self.inputs:
            self.ctx.inputs.settings = self.inputs.settings.get_dict()
        else:
            self.ctx.inputs.settings = {}

        if 'options' in self.inputs:
            self.ctx.inputs._options = self.inputs.options.get_dict()
        else:
            self.ctx.inputs._options = {}

        if 'vdw_table' in self.inputs:
            self.ctx.inputs.vdw_table = self.inputs.vdw_table

        # Either automatic_parallelization or options has to be specified
        if not any([key in self.inputs for key in ['options', 'automatic_parallelization']]):
            self.abort_nowait('you have to specify either the options or automatic_parallelization input')
            return

        # If automatic parallelization is not enabled, we better make sure that the options satisfy minimum requirements
        if 'automatic_parallelization' not in self.inputs:
            num_machines = self.ctx.inputs['_options'].get('resources', {}).get('num_machines', None)
            max_wallclock_seconds = self.ctx.inputs['_options'].get('max_wallclock_seconds', None)

            if num_machines is None or max_wallclock_seconds is None:
                self.abort_nowait("no automatic_parallelization requested, but the options do not specify both '{}' and '{}'"
                    .format('num_machines', 'max_wallclock_seconds'))

        # Validate the inputs related to pseudopotentials
        structure = self.inputs.structure
        pseudos = self.inputs.get('pseudos', None)
        pseudo_family = self.inputs.get('pseudo_family', None)

        try:
            self.ctx.inputs.pseudo = validate_and_prepare_pseudos_inputs(structure, pseudos, pseudo_family)
        except ValueError as exception:
            self.abort_nowait('{}'.format(exception))

    def should_run_init(self):
        """
        Return whether an initialization calculation should be run, which is the case if the user wants
        to use automatic parallelization and has specified the ParameterData node in the inputs
        """
        return 'automatic_parallelization' in self.inputs

    def validate_init_inputs(self):
        """
        Validate the inputs that are required for the initialization calculation. The automatic_parallelization
        input expects a ParameterData node with the following keys:

            * max_wallclock_seconds
            * target_time_seconds
            * max_num_machines

        If any of these keys are not set or any superfluous keys are specified, the workchain will abort.
        """
        automatic_parallelization = self.inputs.automatic_parallelization.get_dict()

        expected_keys = ['max_wallclock_seconds', 'target_time_seconds', 'max_num_machines']
        received_keys = [(key, automatic_parallelization.get(key, None)) for key in expected_keys]
        remaining_keys = [key for key in automatic_parallelization.keys() if key not in expected_keys]

        for k, v in [(key, value) for key, value in received_keys if value is None]:
            self.abort_nowait('required key "{}" in automatic_parallelization input not found'.format(k))
            return

        if remaining_keys:
            self.abort_nowait('detected unrecognized keys in the automatic_parallelization input: {}'.format(' '.join(remaining_keys)))
            return

        # Add the calculation mode to the automatic parallelization dictionary
        self.ctx.automatic_parallelization = {
            'max_wallclock_seconds': automatic_parallelization['max_wallclock_seconds'],
            'target_time_seconds': automatic_parallelization['target_time_seconds'],
            'max_num_machines': automatic_parallelization['max_num_machines'],
            'calculation_mode': self.ctx.inputs.parameters['CONTROL']['calculation']
        }

        self.ctx.inputs._options.setdefault('resources', {})['num_machines'] = automatic_parallelization['max_num_machines']
        self.ctx.inputs._options['max_wallclock_seconds'] = automatic_parallelization['max_wallclock_seconds']

    def run_init(self):
        """
        Run a first dummy pw calculation that will exit straight after initialization. At that
        point it will have generated some general output pertaining to the dimensions of the
        calculation which we can use to distribute available computational resources
        """
        inputs = deepcopy(self.ctx.inputs)

        # Set the initialization flag and the initial default options
        inputs['settings']['ONLY_INITIALIZATION'] = True
        inputs['_options'] = update_mapping(inputs['_options'], get_default_options())

        # Prepare the final input dictionary
        inputs = self._prepare_process_inputs(inputs)
        process = PwCalculation.process()
        running = submit(process, **inputs)

        self.report('launching initialization PwCalculation<{}>'.format(running.pid))

        return ToContext(calculation_init=running)

    def inspect_init(self):
        """
        Use the initialization PwCalculation to determine the required resource settings for the
        requested calculation based on the settings in the automatic_parallelization input
        """
        calculation = self.ctx.calculation_init

        if not calculation.has_finished_ok():
            self.abort_nowait('the initialization calculation did not finish successfully')
            return

        # Get automated parallelization settings
        parallelization = get_pw_parallelization_parameters(calculation, **self.ctx.automatic_parallelization)

        self.report('Determined the following resource settings from automatic_parallelization input: {}'
            .format(parallelization))

        options = self.ctx.inputs._options
        base_resources = options.get('resources', {})
        goal_resources = parallelization['resources']

        scheduler = calculation.get_computer().get_scheduler()
        resources = create_scheduler_resources(scheduler, base_resources, goal_resources)

        cmdline = self.ctx.inputs.settings.get('cmdline', [])
        cmdline = cmdline_remove_npools(cmdline)
        cmdline.extend(['-nk', str(parallelization['npools'])])

        # Set the new cmdline setting and resource options
        self.ctx.inputs.settings['cmdline'] = cmdline
        self.ctx.inputs._options = update_mapping(options, {'resources': resources})

        return

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
        inputs = deepcopy(self.ctx.inputs)

        if self.ctx.restart_calc:
            inputs['parameters']['CONTROL']['restart_mode'] = 'restart'
            inputs['parent_folder'] = self.ctx.restart_calc.out.remote_folder

        # Prepare the final input dictionary
        inputs = self._prepare_process_inputs(inputs)
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
            self.abort_nowait('unexpected state ({}) of PwCalculation<{}>'
                .format(calculation.get_state(), calculation.pk))

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
        Attach the output parameters and retrieved folder of the last calculation to the outputs
        """
        self.report('workchain completed after {} iterations'.format(self.ctx.iteration))
        self.out('output_parameters', self.ctx.restart_calc.out.output_parameters)
        self.out('remote_folder', self.ctx.restart_calc.out.remote_folder)
        self.out('retrieved', self.ctx.restart_calc.out.retrieved)

        if 'output_structure' in self.ctx.restart_calc.out:
            self.out('output_structure', self.ctx.restart_calc.out.output_structure)

        if 'output_band' in self.ctx.restart_calc.out:
            self.out('output_band', self.ctx.restart_calc.out.output_band)

    def clean(self):
        """
        Clean remote folders of the PwCalculations that were run if the clean_workdir parameter was
        set to true in the Workchain inputs
        """
        if not self.inputs.clean_workdir.value:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = []
        calculations = self.ctx.calculations

        try:
            calculations.append(self.ctx.calculation_init)
        except AttributeError:
            pass

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

    def _prepare_process_inputs(self, inputs):
        """
        Prepare the inputs dictionary for a PwCalculation process. The 'max_seconds' setting in the 'CONTROL' card
        of the parameters will be set to a fraction of the 'max_wallclock_seconds' that will be given to the job via
        the '_options' dictionary. This will prevent the job from being prematurely terminated by the scheduler without
        getting the chance to exit cleanly. Any remaining bare dictionaries in the inputs dictionary will be wrapped
        in a ParameterData data node except for the '_options' key which should remain a standard dictionary
        """
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