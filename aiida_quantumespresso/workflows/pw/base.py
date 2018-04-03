# -*- coding: utf-8 -*-
from aiida.orm import Code
from aiida.orm.data.base import Bool, Int, Str
from aiida.orm.data.folder import FolderData
from aiida.orm.data.remote import RemoteData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array import ArrayData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.singlefile import SinglefileData
from aiida.orm.utils import CalculationFactory
from aiida.common.extendeddicts import AttributeDict
from aiida.work.run import submit
from aiida.work.workchain import ToContext, if_, while_
from aiida_quantumespresso.common.exceptions import UnexpectedCalculationFailure
from aiida_quantumespresso.common.workchain.utils import ErrorHandlerReport
from aiida_quantumespresso.common.workchain.utils import register_error_handler
from aiida_quantumespresso.common.workchain.base.restart import BaseRestartWorkChain
from aiida_quantumespresso.utils.defaults.calculation import pw as qe_defaults
from aiida_quantumespresso.utils.mapping import update_mapping
from aiida_quantumespresso.utils.pseudopotential import validate_and_prepare_pseudos_inputs
from aiida_quantumespresso.utils.resources import get_default_options
from aiida_quantumespresso.utils.resources import get_pw_parallelization_parameters
from aiida_quantumespresso.utils.resources import cmdline_remove_npools
from aiida_quantumespresso.utils.resources import create_scheduler_resources


PwCalculation = CalculationFactory('quantumespresso.pw')


class PwBaseWorkChain(BaseRestartWorkChain):
    """
    Base workchain to launch a Quantum Espresso pw.x calculation
    """
    _verbose = True
    _calculation_class = PwCalculation
    _error_handler_entry_point = 'aiida_quantumespresso.workflow_error_handlers.pw.base'

    defaults = AttributeDict({
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
        spec.outline(
            cls.setup,
            cls.validate_inputs,
            if_(cls.should_run_init)(
                cls.validate_init_inputs,
                cls.run_init,
                cls.inspect_init,
            ),
            while_(cls.should_run_calculation)(
                cls.prepare_calculation,
                cls.run_calculation,
                cls.inspect_calculation,
            ),
            cls.results,
        )
        spec.output('output_array', valid_type=ArrayData, required=False)
        spec.output('output_band', valid_type=BandsData, required=False)
        spec.output('output_structure', valid_type=StructureData, required=False)
        spec.output('output_parameters', valid_type=ParameterData)
        spec.output('remote_folder', valid_type=RemoteData)
        spec.output('retrieved', valid_type=FolderData)

    def validate_inputs(self):
        """
        Validate inputs that depend might depend on each other and cannot be validated by the spec. Also define
        dictionary `inputs` in the context, that will contain the inputs for the calculation that will be launched
        in the `run_calculation` step.
        """
        self.ctx.inputs = AttributeDict({
            'code': self.inputs.code,
            'structure': self.inputs.structure,
            'kpoints': self.inputs.kpoints,
            'parameters': self.inputs.parameters.get_dict()
        })

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
        parallelization = self.inputs.automatic_parallelization.get_dict()

        expected_keys = ['max_wallclock_seconds', 'target_time_seconds', 'max_num_machines']
        received_keys = [(key, parallelization.get(key, None)) for key in expected_keys]
        remaining_keys = [key for key in parallelization.keys() if key not in expected_keys]

        for k, v in [(key, value) for key, value in received_keys if value is None]:
            self.abort_nowait('required key "{}" in automatic_parallelization input not found'.format(k))
            return

        if remaining_keys:
            self.abort_nowait('detected unrecognized keys in the automatic_parallelization input: {}'
                .format(' '.join(remaining_keys)))
            return

        # Add the calculation mode to the automatic parallelization dictionary
        self.ctx.automatic_parallelization = {
            'max_wallclock_seconds': parallelization['max_wallclock_seconds'],
            'target_time_seconds': parallelization['target_time_seconds'],
            'max_num_machines': parallelization['max_num_machines'],
            'calculation_mode': self.ctx.inputs.parameters['CONTROL']['calculation']
        }

        self.ctx.inputs._options.setdefault('resources', {})['num_machines'] = parallelization['max_num_machines']
        self.ctx.inputs._options['max_wallclock_seconds'] = parallelization['max_wallclock_seconds']

    def run_init(self):
        """
        Run a first dummy pw calculation that will exit straight after initialization. At that
        point it will have generated some general output pertaining to the dimensions of the
        calculation which we can use to distribute available computational resources
        """
        inputs = self.ctx.inputs

        # Set the initialization flag and the initial default options
        inputs.settings['ONLY_INITIALIZATION'] = True
        inputs._options = update_mapping(inputs['_options'], get_default_options())

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

        node = ParameterData(dict=parallelization)
        self.out('automatic_parallelization', node)
        self.report('results of automatic parallelization in {}<{}>'.format(node.__class__.__name__, node.pk))

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

        # Remove the only initialization flag
        self.ctx.inputs.settings.pop('ONLY_INITIALIZATION')

        return

    def prepare_calculation(self):
        """
        Prepare the inputs for the next calculation
        """
        if self.ctx.restart_calc:
            self.ctx.inputs.parameters['CONTROL']['restart_mode'] = 'restart'
            self.ctx.inputs.parent_folder = self.ctx.restart_calc.out.remote_folder

    def _prepare_process_inputs(self, inputs):
        """
        The 'max_seconds' setting in the 'CONTROL' card of the parameters will be set to a fraction of the
        'max_wallclock_seconds' that will be given to the job via the '_options' dictionary. This will prevent the job
        from being prematurely terminated by the scheduler without getting the chance to exit cleanly.
        """
        max_wallclock_seconds = inputs._options['max_wallclock_seconds']
        max_seconds_factor = self.defaults.delta_factor_max_seconds
        max_seconds = max_wallclock_seconds * max_seconds_factor
        inputs.parameters['CONTROL']['max_seconds'] = max_seconds

        return super(PwBaseWorkChain, self)._prepare_process_inputs(inputs)


@register_error_handler(PwBaseWorkChain, 500)
def _handle_error_read_namelists(self, calculation):
    """
    The calculation failed because it could not read the generated input file
    """
    if any(['read_namelists' in w for w in calculation.res.warnings]):
        self.abort_nowait('PwCalculation<{}> failed because of an invalid input file'.format(calculation.pk))
        return ErrorHandlerReport(True, True)


@register_error_handler(PwBaseWorkChain, 400)
def _handle_error_exceeded_maximum_walltime(self, calculation):
    """
    Calculation ended nominally but ran out of allotted wall time
    """
    if 'Maximum CPU time exceeded' in calculation.res.warnings:
        self.ctx.restart_calc = calculation
        self.report('PwCalculation<{}> terminated because maximum wall time was exceeded, restarting'
            .format(calculation.pk))
        return ErrorHandlerReport(True, True)


@register_error_handler(PwBaseWorkChain, 300)
def _handle_error_diagonalization(self, calculation):
    """
    Diagonalization failed with current scheme. Try to restart from previous clean calculation with different scheme
    """
    input_parameters = calculation.inp.parameters.get_dict()
    input_electrons = input_parameters.get('ELECTRONS', {})
    diagonalization = input_electrons.get('diagonalization', self.defaults['qe']['diagonalization'])

    if ((
        any(['too many bands are not converged' in w for w in calculation.res.warnings]) or
        any(['eigenvalues not converged' in w for w in calculation.res.warnings])
    ) and (
        diagonalization == 'david'
    )):
        new_diagonalization = 'cg'
        self.ctx.inputs.parameters['ELECTRONS']['diagonalization'] = 'cg'
        self.ctx.restart_calc = calculation
        self.report('PwCalculation<{}> failed to diagonalize with "{}" scheme'.format(diagonalization))
        self.report('Restarting with diagonalization scheme "{}"'.format(new_diagonalization))
        return ErrorHandlerReport(True, True)


@register_error_handler(PwBaseWorkChain, 200)
def _handle_error_convergence_not_reached(self, calculation):
    """
    At the end of the scf cycle, the convergence threshold was not reached. We simply restart
    from the previous calculation without changing any of the input parameters
    """
    if 'The scf cycle did not reach convergence.' in calculation.res.warnings:
        self.ctx.restart_calc = calculation
        self.report('PwCalculation<{}> did not converge, restart from previous calculation'.format(calculation.pk))
        return ErrorHandlerReport(True, True)


@register_error_handler(PwBaseWorkChain, 100)
def _handle_error_unrecognized_by_parser(self, calculation):
    """
    Calculation failed with an error that was not recognized by the parser and was attached
    wholesale to the warnings. We treat it as an unexpected failure and raise the exception
    """
    warnings = calculation.res.warnings
    if (any(['%%%' in w for w in warnings]) or any(['Error' in w for w in warnings])):
        raise UnexpectedCalculationFailure('PwCalculation<{}> failed due to an unknown reason'.format(calculation.pk))
