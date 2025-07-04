# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO xspectra.x calculation with automated error handling and restarts."""
import warnings

from aiida import orm
from aiida.common import AttributeDict
from aiida.common.lang import type_check
from aiida.engine import BaseRestartWorkChain, ProcessHandlerReport, process_handler, while_
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import create_kpoints_from_distance
from aiida_quantumespresso.workflows.protocols.utils import ProtocolMixin

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    XspectraCalculation = CalculationFactory('quantumespresso.xspectra')

warnings.warn(
    'This module is deprecated and will be removed soon as part of migrating XAS and XPS workflows to a new repository.'
    '\nThe new repository can be found at: https://github.com/aiidaplugins/aiida-qe-xspec.', FutureWarning
)


class XspectraBaseWorkChain(ProtocolMixin, BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO xspectra.x calculation with automated error handling and restarts."""

    _process_class = XspectraCalculation

    defaults = AttributeDict({'delta_factor_time_limit': 0.90})

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        super().define(spec)
        spec.expose_inputs(XspectraCalculation, namespace='xspectra', exclude=('kpoints',))
        spec.input(
            'kpoints',
            valid_type=orm.KpointsData,
            required=False,
            help='An explicit k-points mesh. Either this or `kpoints_distance` has to be provided.'
        )
        spec.input(
            'kpoints_distance',
            valid_type=orm.Float,
            required=False,
            help='The minimum desired distance in 1/Å between k-points in reciprocal space. The explicit k-points will '
            'be generated automatically by a calculation function based on the input structure.'
        )
        spec.input(
            'kpoints_force_parity',
            valid_type=orm.Bool,
            required=False,
            help='Optional input when constructing the k-points based on a desired `kpoints_distance`. Setting this to '
            '`True` will force the k-point mesh to have an even number of points along each lattice vector except '
            'for any non-periodic directions.'
        )
        spec.expose_outputs(XspectraCalculation)
        spec.outline(
            cls.setup,
            cls.validate_kpoints,
            while_(cls.should_run_process)(
                cls.prepare_process,
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )
        spec.exit_code(
            202,
            'ERROR_INVALID_INPUT_KPOINTS',
            message='Neither the `kpoints` nor the `kpoints_distance` input was specified.'
        )
        spec.exit_code(
            300, 'ERROR_UNRECOVERABLE_FAILURE', message='The calculation failed with an unrecoverable error.'
        )

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        from ..protocols import xspectra as xs_protocols
        return files(xs_protocols) / 'base.yaml'

    @classmethod
    def get_builder_from_protocol(
        cls, code, core_wfc_data, parent_folder, abs_atom_marker='X', protocol=None, overrides=None, options=None, **_
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param code: the ``Code`` instance configured for the ``quantumespresso.xspectra`` plugin.
        :param core_wfc_data: a ``SinglefileData`` object for the initial-state core
                              wavefunction, normally derived from upf2plotcore.sh, required
                              for the xspectra.x calculation.
        :param parent_folder: a ``RemoteData`` object for the parent calculation (either pw.x
                              or xspectra.x).
        :param abs_atom_marker: the name given to the Kind representing the absorbing atom.
                                Matches to Kind.name for the absorbing atom in the structure,
                                used to set the parameter `xiabs` automatically.
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

        metadata = inputs['xspectra']['metadata']
        parameters = inputs['xspectra']['parameters']

        if options:
            metadata['options'] = recursive_merge(inputs['xspectra']['metadata']['options'], options)

        # pylint: disable=no-member
        builder = cls.get_builder()
        parent_calc = parent_folder.creator
        if parent_calc.process_type == 'aiida.calculations:quantumespresso.xspectra':
            builder.kpoints = parent_calc.inputs.kpoints
        elif parent_calc.process_type == 'aiida.calculations:quantumespresso.pw':
            structure = parent_calc.inputs.structure
            kinds_present = sorted([kind.name for kind in structure.kinds])
            if 'kpoints' in inputs['xspectra']:
                builder.kpoints = inputs['xspectra']['kpoints']
            else:
                builder.kpoints_distance = orm.Float(inputs['kpoints_distance'])
            builder.kpoints_force_parity = orm.Bool(inputs['kpoints_force_parity'])
            if abs_atom_marker not in kinds_present:
                raise ValueError(f'given abs_atom_marker `{abs_atom_marker}` does not match a kind in the structure.')
            else:
                parameters['INPUT_XSPECTRA']['xiabs'] = kinds_present.index(abs_atom_marker) + 1
        else:
            raise ValueError(f'process type of `parent_folder` creator `{parent_calc.process_type}` is not supported.')

        builder.xspectra['code'] = code
        builder.xspectra['parent_folder'] = parent_folder
        builder.xspectra['core_wfc_data'] = core_wfc_data
        builder.xspectra['parameters'] = orm.Dict(parameters)
        builder.xspectra['metadata'] = metadata
        if 'settings' in inputs['xspectra']:
            builder.xspectra['settings'] = orm.Dict(inputs['xspectra']['settings'])
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])
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
        self.ctx.inputs = AttributeDict(self.exposed_inputs(XspectraCalculation, 'xspectra'))

        self.ctx.inputs.parameters = self.ctx.inputs.parameters.get_dict()

    def validate_kpoints(self):
        """Validate the inputs related to k-points.

        Either an explicit `KpointsData` with given mesh, or a desired k-points distance should be specified. In
        the case of the latter, the `KpointsData` will be constructed for the input `StructureData` using the
        `create_kpoints_from_distance` calculation function.
        """
        if all(key not in self.inputs for key in ['kpoints', 'kpoints_distance']):
            return self.exit_codes.ERROR_INVALID_INPUT_KPOINTS

        # xspectra.x can only work with a kpoints mesh, so we check that the KpointsData has a
        # mesh property:
        if 'kpoints' in self.inputs:
            try:
                self.inputs.kpoints.get_kpoints_mesh()
            except AttributeError as exception:
                raise AttributeError('XSpectra calculations cannot use an explicit kpoints list.') from exception

        try:
            kpoints = self.inputs.kpoints
        except AttributeError:
            inputs = {
                'structure': self.inputs.xspectra.parent_folder.creator.inputs.structure,
                'distance': self.inputs.kpoints_distance,
                'force_parity': self.inputs.get('kpoints_force_parity', orm.Bool(False)),
                'metadata': {
                    'call_link_label': 'create_kpoints_from_distance'
                }
            }
            kpoints = create_kpoints_from_distance(**inputs)  # pylint: disable=unexpected-keyword-arg

        self.ctx.inputs.kpoints = kpoints

    def set_max_seconds(self, max_wallclock_seconds):
        """Set the `max_seconds` to a fraction of `max_wallclock_seconds` option to prevent out-of-walltime problems.

        :param max_wallclock_seconds: the maximum wallclock time that will be set in the scheduler settings.
        """
        max_seconds_factor = self.defaults.delta_factor_time_limit
        max_seconds = max_wallclock_seconds * max_seconds_factor
        self.ctx.inputs.parameters['INPUT_XSPECTRA']['time_limit'] = max_seconds

    def prepare_process(self):
        """Prepare the inputs for the next calculation.

        If a `restart_calc` has been set in the context, its `remote_folder` will be used as the `parent_folder` input
        for the next calculation and the `restart_mode` is set to `restart`.
        """
        max_wallclock_seconds = self.ctx.inputs.metadata.options.get('max_wallclock_seconds', None)

        if max_wallclock_seconds is not None and 'time_limit' not in self.ctx.inputs.parameters['INPUT_XSPECTRA']:
            self.set_max_seconds(max_wallclock_seconds)

        if self.ctx.restart_calc:
            self.ctx.inputs.parameters['INPUT_XSPECTRA']['restart_mode'] = 'restart'
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

    @process_handler(priority=610, exit_codes=XspectraCalculation.exit_codes.ERROR_SCHEDULER_OUT_OF_WALLTIME)
    def handle_scheduler_out_of_walltime(self, node):
        """Handle `ERROR_SCHEDULER_OUT_OF_WALLTIME` exit code: decrease the time_limit and restart from scratch."""

        # Decrease `max_seconds` significantly in order to make sure that the calculation has the time to shut down
        # neatly before reaching the scheduler wall time and one can restart from this calculation.
        factor = 0.5
        max_seconds = self.ctx.inputs.parameters.get('INPUT_XSPECTRA', {}).get('time_limit', None)
        if max_seconds is None:
            max_seconds = self.ctx.inputs.metadata.options.get(
                'max_wallclock_seconds', None
            ) * self.defaults.delta_factor_time_limit
        max_seconds_new = max_seconds * factor

        self.ctx.restart_calc = node
        self.ctx.inputs.parameters.setdefault('INPUT_XSPECTRA', {})['restart_type'] = 'from_scratch'
        self.ctx.inputs.parameters.setdefault('INPUT_XSPECTRA', {})['time_limit'] = max_seconds_new

        action = f'reduced max_seconds from {max_seconds} to {max_seconds_new} and restarting'
        self.report_error_handled(node, action)
        return ProcessHandlerReport(True)

    @process_handler(priority=580, exit_codes=XspectraCalculation.exit_codes.ERROR_OUT_OF_WALLTIME)
    def handle_out_of_walltime(self, node):
        """Handle `ERROR_OUT_OF_WALLTIME` exit code: calculation shut down neatly and we can simply restart."""
        self.ctx.restart_calc = node
        self.report_error_handled(node, 'simply restart from the last calculation')
        return ProcessHandlerReport(True)
