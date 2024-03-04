# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO ph.x calculation with automated error handling and restarts."""
from typing import Mapping

from aiida import orm
from aiida.common import AttributeDict
from aiida.common.lang import type_check
from aiida.engine import BaseRestartWorkChain, ProcessHandlerReport, process_handler, while_
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import create_kpoints_from_distance
from aiida_quantumespresso.calculations.functions.merge_ph_outputs import merge_ph_outputs
from aiida_quantumespresso.common.types import ElectronicType
from aiida_quantumespresso.workflows.protocols.utils import ProtocolMixin

PhCalculation = CalculationFactory('quantumespresso.ph')
PwCalculation = CalculationFactory('quantumespresso.pw')


class PhBaseWorkChain(ProtocolMixin, BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO ph.x calculation with automated error handling and restarts."""

    _process_class = PhCalculation

    defaults = AttributeDict({
        'delta_factor_max_seconds': 0.95,
        'delta_factor_alpha_mix': 0.90,
        'alpha_mix': 0.70,
    })

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.expose_inputs(PhCalculation, namespace='ph', exclude=('qpoints', ))
        spec.input('only_initialization', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.input('qpoints', valid_type=orm.KpointsData, required=False,
            help='An explicit qpoints list or mesh. Either this or `qpoints_distance` should to be provided.')
        spec.input('qpoints_distance', valid_type=orm.Float, required=False,
            help='The minimum desired distance in 1/Å between qpoints in reciprocal space. The explicit qpoints will '
                 'be generated automatically by a calculation function based on the input structure.')
        spec.input('qpoints_force_parity', valid_type=orm.Bool, required=False,
            help='Optional input when constructing the qpoints based on a desired `qpoints_distance`. Setting this to '
                 '`True` will force the qpoint mesh to have an even number of points along each lattice vector except '
                 'for any non-periodic directions.')
        spec.inputs.validator = cls.validate_inputs
        spec.outline(
            cls.setup,
            cls.validate_parameters,
            cls.set_qpoints,
            while_(cls.should_run_process)(
                cls.prepare_process,
                cls.run_process,
                cls.inspect_process,
            ),
            cls.create_merged_output,
            cls.results,
        )
        spec.expose_outputs(PhCalculation, exclude=('retrieved_folder',))
        spec.exit_code(204, 'ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED',
            message='The `metadata.options` did not specify both `resources.num_machines` and `max_wallclock_seconds`. '
                    'This exit status has been deprecated as the check it corresponded to was incorrect.')
        spec.exit_code(300, 'ERROR_UNRECOVERABLE_FAILURE',
            message='The calculation failed with an unrecoverable error.')
        spec.exit_code(401, 'ERROR_MERGING_QPOINTS',
            message='The work chain failed to merge the q-points data from multiple `PhCalculation`s because not all '
                    'q-points were parsed.')
        # yapf: enable

    @classmethod
    def validate_inputs(cls, value, port_namespace):  # pylint: disable=unused-argument
        """Validate the top level namespace."""

        if (('qpoints_distance' in port_namespace or 'qpoints' in port_namespace) and
            'qpoints_distance' not in value and 'qpoints' not in value):
            return 'Neither `qpoints` nor `qpoints_distance` were specified.'

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        from ..protocols import ph as ph_protocols
        return files(ph_protocols) / 'base.yaml'

    @classmethod
    def get_builder_from_protocol(
        cls,
        code,
        parent_folder=None,
        protocol=None,
        overrides=None,
        electronic_type=ElectronicType.METAL,
        options=None,
        **_
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param code: the ``Code`` instance configured for the ``quantumespresso.ph`` plugin.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param options: A dictionary of options that will be recursively set for the ``metadata.options`` input of all
            the ``CalcJobs`` that are nested in this work chain.
        :param electronic_type: indicate the electronic character of the system through ``ElectronicType`` instance.
        :return: a process builder instance with all inputs defined ready for launch.
        """
        from aiida_quantumespresso.workflows.protocols.utils import recursive_merge

        if isinstance(code, str):
            code = orm.load_code(code)

        type_check(code, orm.AbstractCode)
        type_check(electronic_type, ElectronicType)

        if electronic_type not in [ElectronicType.METAL, ElectronicType.INSULATOR]:
            raise NotImplementedError(f'electronic type `{electronic_type}` is not supported.')

        inputs = cls.get_protocol_inputs(protocol, overrides)

        if electronic_type is ElectronicType.INSULATOR:
            inputs['ph']['parameters']['INPUTPH']['epsil'] = True

        metadata = inputs['ph']['metadata']

        if options:
            metadata['options'] = recursive_merge(inputs['ph']['metadata']['options'], options)

        # pylint: disable=no-member
        builder = cls.get_builder()
        builder.ph['code'] = code
        if parent_folder is not None:
            builder.ph['parent_folder'] = parent_folder
        builder.ph['parameters'] = orm.Dict(inputs['ph']['parameters'])
        builder.ph['metadata'] = metadata
        if 'settings' in inputs['ph']:
            builder.ph['settings'] = orm.Dict(inputs['ph']['settings'])
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])

        if 'qpoints' in inputs:
            qpoints_mesh = inputs['qpoints']
            qpoints = orm.KpointsData()
            qpoints.set_kpoints_mesh(qpoints_mesh)
            builder.qpoints = qpoints
        else:
            builder.qpoints_distance = orm.Float(inputs['qpoints_distance'])
            builder.qpoints_force_parity = orm.Bool(inputs['qpoints_force_parity'])

        builder.max_iterations = orm.Int(inputs['max_iterations'])
        # pylint: enable=no-member

        return builder

    def setup(self):
        """Call the `setup` of the `BaseRestartWorkChain` and then create the inputs dictionary in `self.ctx.inputs`.

        This `self.ctx.inputs` dictionary will be used by the `BaseRestartWorkChain` to submit the calculations in the
        internal loop.
        """
        super().setup()
        self.ctx.restart_calc = None
        self.ctx.inputs = AttributeDict(self.exposed_inputs(PhCalculation, 'ph'))

    def validate_parameters(self):
        """Validate inputs that might depend on each other and cannot be validated by the spec."""
        self.ctx.inputs.parameters = self.ctx.inputs.parameters.get_dict()
        self.ctx.inputs.settings = self.ctx.inputs.settings.get_dict() if 'settings' in self.ctx.inputs else {}

        self.ctx.inputs.parameters.setdefault('INPUTPH', {})
        self.ctx.inputs.parameters['INPUTPH']['recover'] = 'parent_folder' in self.ctx.inputs

        if self.inputs.only_initialization.value:
            self.ctx.inputs.settings['ONLY_INITIALIZATION'] = True

    def set_qpoints(self):
        """Set the inputs related to qpoints.

        Either an explicit `KpointsData` with given mesh/path, or a desired qpoints distance should be specified. In
        the case of the latter, the `KpointsData` will be constructed for the input `StructureData`
        from the parent_folder using the `create_kpoints_from_distance` calculation function.
        """
        try:
            qpoints = self.inputs.qpoints
        except AttributeError:

            try:
                structure = self.ctx.inputs.parent_folder.creator.output.output_structure
            except AttributeError:
                structure = self.ctx.inputs.parent_folder.creator.inputs.structure

            inputs = {
                'structure': structure,
                'distance': self.inputs.qpoints_distance,
                'force_parity': self.inputs.get('qpoints_force_parity', orm.Bool(False)),
                'metadata': {
                    'call_link_label': 'create_qpoints_from_distance'
                }
            }
            qpoints = create_kpoints_from_distance(**inputs)

        self.ctx.inputs['qpoints'] = qpoints

    def set_max_seconds(self, max_wallclock_seconds: None):
        """Set the `max_seconds` to a fraction of `max_wallclock_seconds` option to prevent out-of-walltime problems.

        :param max_wallclock_seconds: the maximum wallclock time that will be set in the scheduler settings.
        """
        max_seconds_factor = self.defaults.delta_factor_max_seconds
        max_seconds = max_wallclock_seconds * max_seconds_factor
        self.ctx.inputs.parameters['INPUTPH']['max_seconds'] = max_seconds

    def prepare_process(self):
        """Prepare the inputs for the next calculation.

        If a `restart_calc` has been set in the context, its `remote_folder` will be used as the `parent_folder` input
        for the next calculation and the `restart_mode` is set to `restart`.
        """
        max_wallclock_seconds = self.ctx.inputs.metadata.options.get('max_wallclock_seconds', None)

        if max_wallclock_seconds is not None and 'max_seconds' not in self.ctx.inputs.parameters['INPUTPH']:
            self.set_max_seconds(max_wallclock_seconds)

        if self.ctx.restart_calc:
            self.ctx.inputs.parameters['INPUTPH']['recover'] = True
            self.ctx.inputs.parent_folder = self.ctx.restart_calc.outputs.remote_folder

    def create_merged_output(self):
        """Merge outputs from multiple ``PhCalculation`` runs called by the workchain if necessary."""
        if self.inputs.only_initialization.value:
            return

        if len(self.ctx.children) == 1:
            return

        output_dict = {
            'output_' + str(index + 1): child.outputs.output_parameters
            for index, child in enumerate(self.ctx.children)
        }

        num_qpoints = self.ctx.children[0].outputs.output_parameters['number_of_qpoints']
        num_qpoints = self.ctx.inputs.parameters['INPUTPH'].get('last_q', num_qpoints) \
            - self.ctx.inputs.parameters['INPUTPH'].get('start_q', 1) + 1
        num_qpoints_found = sum(
            len(output['number_of_irr_representations_for_each_q']) for output in output_dict.values()
        )

        if num_qpoints_found == num_qpoints:
            self.report(f'Merging {num_qpoints} q-points data from `PhCalculation`s.')
            self.ctx.merged_output_parameters = merge_ph_outputs(**output_dict)
        else:
            self.report(f'Only {num_qpoints_found} of {num_qpoints} q-points were parsed.')
            return self.exit_codes.ERROR_MERGING_QPOINTS

    def get_outputs(self, node) -> Mapping[str, orm.Node]:
        """Return a mapping of the outputs that should be attached as outputs to the work chain."""
        outputs = super().get_outputs(node)

        if 'merged_output_parameters' in self.ctx:
            outputs['output_parameters'] = self.ctx.merged_output_parameters

        return outputs

    def report_error_handled(self, calculation, action):
        """Report an action taken for a calculation that has failed.

        This should be called in a registered error handler if its condition is met and an action was taken.

        :param calculation: the failed calculation node
        :param action: a string message with the action taken
        """
        arguments = [calculation.process_label, calculation.pk, calculation.exit_status, calculation.exit_message]
        self.report('{}<{}> failed with exit status {}: {}'.format(*arguments))
        self.report(f'Action taken: {action}')

    @process_handler(priority=600)
    def handle_unrecoverable_failure(self, node):
        """Handle calculations with an exit status below 400 which are unrecoverable, so abort the work chain."""
        if node.is_failed and node.exit_status < 400:
            self.report_error_handled(node, 'unrecoverable error, aborting...')
            return ProcessHandlerReport(True, self.exit_codes.ERROR_UNRECOVERABLE_FAILURE)

    @process_handler(priority=610, exit_codes=PhCalculation.exit_codes.ERROR_SCHEDULER_OUT_OF_WALLTIME)
    def handle_scheduler_out_of_walltime(self, node):
        """Handle `ERROR_SCHEDULER_OUT_OF_WALLTIME` exit code: decrease the max_secondes and restart from scratch."""

        # Decrease `max_seconds` significantly in order to make sure that the calculation has the time to shut down
        # neatly before reaching the scheduler wall time and one can restart from this calculation.
        factor = 0.5
        max_seconds = self.ctx.inputs.parameters.get('INPUTPH', {}).get('max_seconds', None)
        if max_seconds is None:
            max_seconds = self.ctx.inputs.metadata.options.get(
                'max_wallclock_seconds', None
            ) * self.defaults.delta_factor_max_seconds
        max_seconds_new = max_seconds * factor

        self.ctx.restart_calc = node
        self.ctx.inputs.parameters.setdefault('INPUTPH', {})['recover'] = False
        self.ctx.inputs.parameters.setdefault('INPUTPH', {})['max_seconds'] = max_seconds_new

        action = f'reduced max_seconds from {max_seconds} to {max_seconds_new} and restarting'
        self.report_error_handled(node, action)
        return ProcessHandlerReport(True)

    @process_handler(priority=585, exit_codes=PhCalculation.exit_codes.ERROR_COMPUTING_CHOLESKY)
    def handle_diagonalization_errors(self, calculation):
        """Handle known issues related to the diagonalization.

        Switch to ``diagonalization = 'cg'`` if not already running with this setting, and restart from the charge
        density. In case the run already used conjugate gradient diagonalization, abort.
        """
        if self.ctx.inputs.parameters['INPUTPH'].get('diagonalization', None) == 'cg':
            action = 'found diagonalization issues but already running with conjugate gradient algorithm, aborting...'
            self.report_error_handled(calculation, action)
            return ProcessHandlerReport(True, self.exit_codes.ERROR_UNRECOVERABLE_FAILURE)

        self.ctx.inputs.parameters['INPUTPH']['diagonalization'] = 'cg'
        action = 'found diagonalization issues, switching to conjugate gradient diagonalization.'
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True)

    @process_handler(priority=580, exit_codes=PhCalculation.exit_codes.ERROR_OUT_OF_WALLTIME)
    def handle_out_of_walltime(self, node):
        """Handle `ERROR_OUT_OF_WALLTIME` exit code: calculation shut down neatly and we can simply restart."""
        self.ctx.restart_calc = node
        self.report_error_handled(node, 'simply restart from the last calculation')
        return ProcessHandlerReport(True)

    @process_handler(priority=410, exit_codes=PhCalculation.exit_codes.ERROR_CONVERGENCE_NOT_REACHED)
    def handle_convergence_not_reached(self, node):
        """Handle `ERROR_CONVERGENCE_NOT_REACHED` exit code: decrease the mixing beta and restart."""
        factor = self.defaults.delta_factor_alpha_mix
        alpha_mix = self.ctx.inputs.parameters.get('INPUTPH', {}).get('alpha_mix(1)', self.defaults.alpha_mix)
        alpha_mix_new = alpha_mix * factor

        self.ctx.restart_calc = node
        self.ctx.inputs.parameters.setdefault('INPUTPH', {})['alpha_mix(1)'] = alpha_mix_new

        action = f'reduced alpha_mix from {alpha_mix} to {alpha_mix_new} and restarting'
        self.report_error_handled(node, action)
        return ProcessHandlerReport(True)
