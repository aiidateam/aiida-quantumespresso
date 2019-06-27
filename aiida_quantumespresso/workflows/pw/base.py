# -*- coding: utf-8 -*-
from __future__ import absolute_import

from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import ToContext, if_, while_
from aiida.plugins import CalculationFactory
from aiida_quantumespresso.common.exceptions import UnexpectedCalculationFailure
from aiida_quantumespresso.common.workchain.utils import ErrorHandlerReport
from aiida_quantumespresso.common.workchain.utils import register_error_handler
from aiida_quantumespresso.common.workchain.base.restart import BaseRestartWorkChain
from aiida_quantumespresso.utils.defaults.calculation import pw as qe_defaults
from aiida_quantumespresso.utils.mapping import update_mapping, prepare_process_inputs
from aiida_quantumespresso.utils.pseudopotential import validate_and_prepare_pseudos_inputs
from aiida_quantumespresso.utils.resources import get_default_options
from aiida_quantumespresso.utils.resources import get_pw_parallelization_parameters
from aiida_quantumespresso.utils.resources import cmdline_remove_npools
from aiida_quantumespresso.utils.resources import create_scheduler_resources
from aiida_quantumespresso.workflows.functions.create_kpoints_from_distance import create_kpoints_from_distance


PwCalculation = CalculationFactory('quantumespresso.pw')


class PwBaseWorkChain(BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO pw.x calculation with automated error handling and restarts"""

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
        # yapf: disable
        super(PwBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=orm.Code,
            help='The code to use for the calculation which is configured to a Quantum ESPRESSO `pw.x` executable.')
        spec.input('structure', valid_type=orm.StructureData,
            help='The input structure.')
        spec.input('kpoints', valid_type=orm.KpointsData, required=False,
            help='An explicit k-points list or mesh. Either this or `kpoints_distance` has to be provided.')
        spec.input('kpoints_distance', valid_type=orm.Float, required=False,
            help='The minimum desired distance in 1/Å between k-points in reciprocal space. The explicit k-points will '
                 'be generated automatically by a calculation function based on the input structure.')
        spec.input('kpoints_force_parity', valid_type=orm.Bool, required=False,
            help='Optional input when constructing the k-points based on a desired `kpoints_distance`. Setting this to '
                 '`True` will force the k-point mesh to have an even number of points along each lattice vector except '
                 'for any non-periodic directions.')
        spec.input('parameters', valid_type=orm.Dict,
            help='The input parameters that are to be used to construct the input file for `pw.x`.')
        spec.input_namespace('pseudos', valid_type=orm.UpfData, required=False, dynamic=True,
            help='A mapping of `UpfData` nodes onto the kind name to which they should apply.')
        spec.input('pseudo_family', valid_type=orm.Str, required=False,
            help='An alternative to specifying the pseudo potentials manually in `pseudos`: one can specify the name '
                 'of an existing pseudo potential family and the work chain will generate the pseudos automatically '
                 'based on the input structure.')
        spec.input('parent_folder', valid_type=orm.RemoteData, required=False,
            help='An optional working directory of a previously completed calculation to restart from.')
        spec.input('vdw_table', valid_type=orm.SinglefileData, required=False,
            help='Optional van der Waals table contained in a `SinglefileData`.')
        spec.input('settings', valid_type=orm.Dict, required=False,
            help='Optional parameters to affect the way the calculation job and the parsing are performed.')
        spec.input('options', valid_type=orm.Dict, required=False,
            help='The metadata options that are to be passed to the calculation job.')
        spec.input('automatic_parallelization', valid_type=orm.Dict, required=False,
            help='When defined, the work chain will first launch an initialization calculation to determine the '
                 'dimensions of the problem, and based on this it will try to set optimal parallelization flags.')

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

        spec.output('automatic_parallelization', valid_type=orm.Dict, required=False,
            help='The results of the automatic parallelization analysis if performed.')
        spec.output('output_array', valid_type=orm.ArrayData, required=False,
            help='The `output_array` output node of the successful calculation if present.')
        spec.output('output_band', valid_type=orm.BandsData, required=False,
            help='The `output_band` output node of the successful calculation if present.')
        spec.output('output_structure', valid_type=orm.StructureData, required=False,
            help='The `output_structure` output node of the successful calculation if present.')
        spec.output('output_parameters', valid_type=orm.Dict,
            help='The `output_parameters` output node of the successful calculation.')
        spec.output('remote_folder', valid_type=orm.RemoteData,
            help='The `remote_folder` output node of the successful calculation.')

        spec.exit_code(301, 'ERROR_INVALID_INPUT_PSEUDO_POTENTIALS',
            message='The explicit `pseudos` or `pseudo_family` could not be used to get the necessary pseudos.')
        spec.exit_code(302, 'ERROR_INVALID_INPUT_KPOINTS',
            message='Neither the `kpoints` nor the `kpoints_distance` input was specified.')
        spec.exit_code(303, 'ERROR_INVALID_INPUT_RESOURCES',
            message='Neither the `options` nor `automatic_parallelization` input was specified.')
        spec.exit_code(304, 'ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED',
            message='The `options` do not specify both `num_machines` and `max_wallclock_seconds`.')
        spec.exit_code(310, 'ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_MISSING_KEY',
            message='Required key for `automatic_parallelization` was not specified.')
        spec.exit_code(311, 'ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_UNRECOGNIZED_KEY',
            message='Unrecognized keys were specified for `automatic_parallelization`.')
        spec.exit_code(401, 'ERROR_INITIALIZATION_CALCULATION_FAILED',
            message='The initialization calculation failed.')
        spec.exit_code(402, 'ERROR_CALCULATION_INVALID_INPUT_FILE',
            message='The calculation failed because it had an invalid input file.')

    def validate_inputs(self):
        """Validate inputs that might depend on each other and cannot be validated by the spec.

        Also define dictionary `inputs` in the context, that will contain the inputs for the calculation that will be
        launched in the `run_calculation` step.
        """
        self.ctx.inputs = AttributeDict({
            'code': self.inputs.code,
            'structure': self.inputs.structure,
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
            self.ctx.inputs.parameters['CONTROL']['restart_mode'] = 'from_scratch'

        if 'settings' in self.inputs:
            self.ctx.inputs.settings = self.inputs.settings.get_dict()
        else:
            self.ctx.inputs.settings = {}

        self.ctx.inputs.metadata = AttributeDict()
        if 'options' in self.inputs:
            self.ctx.inputs.metadata.options = self.inputs.options.get_dict()
        else:
            self.ctx.inputs.metadata.options = {}

        if 'vdw_table' in self.inputs:
            self.ctx.inputs.vdw_table = self.inputs.vdw_table

        # Either automatic_parallelization or options has to be specified
        if 'automatic_parallelization' not in self.inputs and 'options' not in self.ctx.inputs.metadata:
            return self.exit_codes.ERROR_INVALID_INPUT_RESOURCES

        # If automatic parallelization is not enabled, we better make sure that the options satisfy minimum requirements
        if 'automatic_parallelization' not in self.inputs:
            num_machines = self.ctx.inputs.metadata['options'].get('resources', {}).get('num_machines', None)
            max_wallclock_seconds = self.ctx.inputs.metadata['options'].get('max_wallclock_seconds', None)

            if num_machines is None or max_wallclock_seconds is None:
                return self.exit_codes.ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED

            self.set_max_seconds(max_wallclock_seconds)

        # Either a KpointsData with given mesh/path, or a desired distance between k-points should be specified
        if all([key not in self.inputs for key in ['kpoints', 'kpoints_distance']]):
            return self.exit_codes.ERROR_INVALID_INPUT_KPOINTS

        try:
            self.ctx.inputs.kpoints = self.inputs.kpoints
        except AttributeError:
            structure = self.inputs.structure
            distance = self.inputs.kpoints_distance
            force_parity = self.inputs.get('kpoints_force_parity', orm.Bool(False))
            self.ctx.inputs.kpoints = create_kpoints_from_distance(structure, distance, force_parity)

        # Validate the inputs related to pseudopotentials
        structure = self.inputs.structure
        pseudos = self.inputs.get('pseudos', None)
        pseudo_family = self.inputs.get('pseudo_family', None)

        try:
            self.ctx.inputs.pseudos = validate_and_prepare_pseudos_inputs(structure, pseudos, pseudo_family)
        except ValueError as exception:
            self.report('{}'.format(exception))
            return self.exit_codes.ERROR_INVALID_INPUT_PSEUDO_POTENTIALS

    def set_max_seconds(self, max_wallclock_seconds):
        """Set the `max_seconds` to a fraction of `max_wallclock_seconds` option to prevent out-of-walltime problems.

        :param max_wallclock_seconds: the maximum wallclock time that will be set in the scheduler settings.
        """
        max_seconds_factor = self.defaults.delta_factor_max_seconds
        max_seconds = max_wallclock_seconds * max_seconds_factor
        self.ctx.inputs.parameters['CONTROL']['max_seconds'] = max_seconds

    def should_run_init(self):
        """Return whether an initialization calculation should be run.

        :return: boolean, `True` if `automatic_parallelization` was specified in the inputs, `False` otherwise.
        """
        return 'automatic_parallelization' in self.inputs

    def validate_init_inputs(self):
        """Validate the inputs that are required for the initialization calculation.

        The `automatic_parallelization` input expects a `Dict` node with the following keys:

            * max_wallclock_seconds
            * target_time_seconds
            * max_num_machines

        If any of these keys are not set or any superfluous keys are specified, the workchain will abort.
        """
        parallelization = self.inputs.automatic_parallelization.get_dict()

        expected_keys = ['max_wallclock_seconds', 'target_time_seconds', 'max_num_machines']
        received_keys = [(key, parallelization.get(key, None)) for key in expected_keys]
        remaining_keys = [key for key in parallelization.keys() if key not in expected_keys]

        for key, value in [(key, value) for key, value in received_keys if value is None]:
            self.report('required key "{}" in automatic_parallelization input not found'.format(key))
            return self.exit_codes.ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_MISSING_KEY

        if remaining_keys:
            self.report('detected unrecognized keys in the automatic_parallelization input: {}'
                .format(' '.join(remaining_keys)))
            return self.exit_codes.ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_UNRECOGNIZED_KEY

        # Add the calculation mode to the automatic parallelization dictionary
        self.ctx.automatic_parallelization = {
            'max_wallclock_seconds': parallelization['max_wallclock_seconds'],
            'target_time_seconds': parallelization['target_time_seconds'],
            'max_num_machines': parallelization['max_num_machines'],
            'calculation_mode': self.ctx.inputs.parameters['CONTROL']['calculation']
        }

        options = self.ctx.inputs.metadata['options']
        options.setdefault('resources', {})['num_machines'] = parallelization['max_num_machines']
        options['max_wallclock_seconds'] = parallelization['max_wallclock_seconds']

    def run_init(self):
        """Run an initialization `PwCalculation` that will exit after the preamble.

        In the preamble, all the relevant dimensions of the problem are computed which allows us to make an estimate of
        the required resources and what parallelization flags need to be set.
        """
        inputs = self.ctx.inputs

        # Set the initialization flag and the initial default options
        inputs.settings['ONLY_INITIALIZATION'] = True
        inputs.metadata['options'] = update_mapping(inputs.metadata['options'], get_default_options())

        # Prepare the final input dictionary
        inputs = prepare_process_inputs(PwCalculation, inputs)
        running = self.submit(PwCalculation, **inputs)

        self.report('launching initialization PwCalculation<{}>'.format(running.pk))

        return ToContext(calculation_init=running)

    def inspect_init(self):
        """Use the initialization `PwCalculation` to determine the required resource and parallelization settings."""
        calculation = self.ctx.calculation_init

        if not calculation.is_finished_ok:
            return self.exit_codes.ERROR_INITIALIZATION_CALCULATION_FAILED

        # Get automated parallelization settings
        parallelization = get_pw_parallelization_parameters(calculation, **self.ctx.automatic_parallelization)

        # Note: don't do this at home, we are losing provenance here. This should be done by a calculation function
        node = orm.Dict(dict=parallelization).store()
        self.out('automatic_parallelization', node)
        self.report('results of automatic parallelization in {}<{}>'.format(node.__class__.__name__, node.pk))

        options = self.ctx.inputs.metadata['options']
        base_resources = options.get('resources', {})
        goal_resources = parallelization['resources']

        scheduler = calculation.computer.get_scheduler()
        resources = create_scheduler_resources(scheduler, base_resources, goal_resources)

        cmdline = self.ctx.inputs.settings.get('cmdline', [])
        cmdline = cmdline_remove_npools(cmdline)
        cmdline.extend(['-nk', str(parallelization['npools'])])

        # Set the new cmdline setting and resource options
        self.ctx.inputs.settings['cmdline'] = cmdline
        self.ctx.inputs.metadata['options'] = update_mapping(options, {'resources': resources})

        # Remove the only initialization flag
        self.ctx.inputs.settings.pop('ONLY_INITIALIZATION')

        return

    def prepare_calculation(self):
        """Prepare the inputs for the next calculation.

        If a `restart_calc` has been set in the context, its `remote_folder` will be used as the `parent_folder` input
        for the next calculation and the `restart_mode` is set to `restart`.
        """
        if self.ctx.restart_calc:
            self.ctx.inputs.parameters['CONTROL']['restart_mode'] = 'restart'
            self.ctx.inputs.parent_folder = self.ctx.restart_calc.outputs.remote_folder


@register_error_handler(PwBaseWorkChain, 500)
def _handle_error_read_namelists(self, calculation):
    """
    The calculation failed because it could not read the generated input file
    """
    if any(['read_namelists' in w for w in calculation.res.warnings]):
        self.report('PwCalculation<{}> failed because of an invalid input file'.format(calculation.pk))
        return ErrorHandlerReport(True, True, self.exit_codes.ERROR_CALCULATION_INVALID_INPUT_FILE)


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
    input_parameters = calculation.inputs.parameters.get_dict()
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
        self.report('PwCalculation<{}> failed to diagonalize with "{}" scheme'.format(calculation.pk, diagonalization))
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
