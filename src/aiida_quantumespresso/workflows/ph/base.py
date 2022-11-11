# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO ph.x calculation with automated error handling and restarts."""
from aiida import orm
from aiida.common import AttributeDict
from aiida.common.lang import type_check
from aiida.engine import BaseRestartWorkChain, ProcessHandlerReport, process_handler, while_
from aiida.plugins import CalculationFactory

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
        spec.expose_inputs(PhCalculation, namespace='ph')
        spec.input('only_initialization', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.outline(
            cls.setup,
            cls.validate_parameters,
            while_(cls.should_run_process)(
                cls.prepare_process,
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )
        spec.expose_outputs(PhCalculation, exclude=('retrieved_folder',))
        spec.exit_code(204, 'ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED',
            message='The `metadata.options` did not specify both `resources.num_machines` and `max_wallclock_seconds`. '
                    'This exit status has been deprecated as the check it corresponded to was incorrect.')
        spec.exit_code(300, 'ERROR_UNRECOVERABLE_FAILURE',
            message='The calculation failed with an unrecoverable error.')
        # yapf: enable

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        from ..protocols import ph as ph_protocols
        return files(ph_protocols) / 'base.yaml'

    @classmethod
    def get_builder_from_protocol(cls, code, parent_folder=None, protocol=None, overrides=None, options=None, **_):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param code: the ``Code`` instance configured for the ``quantumespresso.ph`` plugin.
        :param structure: the ``StructureData`` instance to use.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param options: A dictionary of options that will be recursively set for the ``metadata.options`` input of all
            the ``CalcJobs`` that are nested in this work chain.
        :return: a process builder instance with all inputs defined ready for launch.
        """
        from aiida_quantumespresso.workflows.protocols.utils import recursive_merge

        if isinstance(code, str):
            code = orm.load_code(code)

        type_check(code, orm.AbstractCode)

        inputs = cls.get_protocol_inputs(protocol, overrides)

        qpoints_mesh = inputs['ph'].pop('qpoints')
        qpoints = orm.KpointsData()
        qpoints.set_kpoints_mesh(qpoints_mesh)

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
        builder.ph['qpoints'] = qpoints
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
