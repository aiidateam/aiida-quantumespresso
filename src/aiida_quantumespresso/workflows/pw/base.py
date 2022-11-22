# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO pw.x calculation with automated error handling and restarts."""
from aiida import orm
from aiida.common import AttributeDict, exceptions
from aiida.common.lang import type_check
from aiida.engine import BaseRestartWorkChain, ExitCode, ProcessHandlerReport, ToContext, if_, process_handler, while_
from aiida.plugins import CalculationFactory, GroupFactory

from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import create_kpoints_from_distance
from aiida_quantumespresso.common.types import ElectronicType, RestartType, SpinType
from aiida_quantumespresso.utils.defaults.calculation import pw as qe_defaults
from aiida_quantumespresso.utils.mapping import prepare_process_inputs, update_mapping
from aiida_quantumespresso.utils.resources import (
    cmdline_remove_npools,
    create_scheduler_resources,
    get_default_options,
    get_pw_parallelization_parameters,
)

from ..protocols.utils import ProtocolMixin

PwCalculation = CalculationFactory('quantumespresso.pw')
SsspFamily = GroupFactory('pseudo.family.sssp')
PseudoDojoFamily = GroupFactory('pseudo.family.pseudo_dojo')
CutoffsPseudoPotentialFamily = GroupFactory('pseudo.family.cutoffs')


class PwBaseWorkChain(ProtocolMixin, BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO pw.x calculation with automated error handling and restarts."""

    # pylint: disable=too-many-public-methods, too-many-statements

    _process_class = PwCalculation

    defaults = AttributeDict({
        'qe': qe_defaults,
        'delta_threshold_degauss': 30,
        'delta_factor_degauss': 0.1,
        'delta_factor_mixing_beta': 0.8,
        'delta_factor_max_seconds': 0.95,
        'delta_factor_nbnd': 0.05,
        'delta_minimum_nbnd': 4,
    })

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.expose_inputs(PwCalculation, namespace='pw', exclude=('kpoints',))
        spec.input('pw.metadata.options.resources', valid_type=dict, required=False)
        spec.input('kpoints', valid_type=orm.KpointsData, required=False,
            help='An explicit k-points list or mesh. Either this or `kpoints_distance` has to be provided.')
        spec.input('kpoints_distance', valid_type=orm.Float, required=False,
            help='The minimum desired distance in 1/Å between k-points in reciprocal space. The explicit k-points will '
                 'be generated automatically by a calculation function based on the input structure.')
        spec.input('kpoints_force_parity', valid_type=orm.Bool, required=False,
            help='Optional input when constructing the k-points based on a desired `kpoints_distance`. Setting this to '
                 '`True` will force the k-point mesh to have an even number of points along each lattice vector except '
                 'for any non-periodic directions.')
        spec.input('automatic_parallelization', valid_type=orm.Dict, required=False,
            help='When defined, the work chain will first launch an initialization calculation to determine the '
                 'dimensions of the problem, and based on this it will try to set optimal parallelization flags.')

        spec.outline(
            cls.setup,
            cls.validate_kpoints,
            if_(cls.should_run_init)(
                cls.validate_init_inputs,
                cls.run_init,
                cls.inspect_init,
            ),
            while_(cls.should_run_process)(
                cls.prepare_process,
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )

        spec.expose_outputs(PwCalculation)
        spec.output('automatic_parallelization', valid_type=orm.Dict, required=False,
            help='The results of the automatic parallelization analysis if performed.')

        spec.exit_code(201, 'ERROR_INVALID_INPUT_PSEUDO_POTENTIALS',
            message='The explicit `pseudos` or `pseudo_family` could not be used to get the necessary pseudos.')
        spec.exit_code(202, 'ERROR_INVALID_INPUT_KPOINTS',
            message='Neither the `kpoints` nor the `kpoints_distance` input was specified.')
        spec.exit_code(203, 'ERROR_INVALID_INPUT_RESOURCES',
            message='Neither the `options` nor `automatic_parallelization` input was specified. '
                    'This exit status has been deprecated as the check it corresponded to was incorrect.')
        spec.exit_code(204, 'ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED',
            message='The `metadata.options` did not specify both `resources.num_machines` and `max_wallclock_seconds`. '
                    'This exit status has been deprecated as the check it corresponded to was incorrect.')
        spec.exit_code(210, 'ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_MISSING_KEY',
            message='Required key for `automatic_parallelization` was not specified.')
        spec.exit_code(211, 'ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_UNRECOGNIZED_KEY',
            message='Unrecognized keys were specified for `automatic_parallelization`.')
        spec.exit_code(300, 'ERROR_UNRECOVERABLE_FAILURE',
            message='The calculation failed with an unidentified unrecoverable error.')
        spec.exit_code(310, 'ERROR_KNOWN_UNRECOVERABLE_FAILURE',
            message='The calculation failed with a known unrecoverable error.')
        spec.exit_code(320, 'ERROR_INITIALIZATION_CALCULATION_FAILED',
            message='The initialization calculation failed.')
        spec.exit_code(501, 'ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF',
            message='Then ionic minimization cycle converged but the thresholds are exceeded in the final SCF.')
        spec.exit_code(710, 'WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED',
            message='The electronic minimization cycle did not reach self-consistency, but `scf_must_converge` '
                    'is `False` and/or `electron_maxstep` is 0.')
        # yapf: enable

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        from ..protocols import pw as pw_protocols
        return files(pw_protocols) / 'base.yaml'

    @classmethod
    def get_builder_from_protocol(
        cls,
        code,
        structure,
        protocol=None,
        overrides=None,
        electronic_type=ElectronicType.METAL,
        spin_type=SpinType.NONE,
        initial_magnetic_moments=None,
        options=None,
        **_
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param structure: the ``StructureData`` instance to use.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param electronic_type: indicate the electronic character of the system through ``ElectronicType`` instance.
        :param spin_type: indicate the spin polarization type to use through a ``SpinType`` instance.
        :param initial_magnetic_moments: optional dictionary that maps the initial magnetic moment of each kind to a
            desired value for a spin polarized calculation. Note that in case the ``starting_magnetization`` is also
            provided in the ``overrides``, this takes precedence over the values provided here. In case neither is
            provided and ``spin_type == SpinType.COLLINEAR``, an initial guess for the magnetic moments is used.
        :param options: A dictionary of options that will be recursively set for the ``metadata.options`` input of all
            the ``CalcJobs`` that are nested in this work chain.
        :return: a process builder instance with all inputs defined ready for launch.
        """
        from aiida_quantumespresso.workflows.protocols.utils import get_starting_magnetization, recursive_merge

        if isinstance(code, str):
            code = orm.load_code(code)

        type_check(code, orm.AbstractCode)
        type_check(electronic_type, ElectronicType)
        type_check(spin_type, SpinType)

        if electronic_type not in [ElectronicType.METAL, ElectronicType.INSULATOR]:
            raise NotImplementedError(f'electronic type `{electronic_type}` is not supported.')

        if spin_type not in [SpinType.NONE, SpinType.COLLINEAR]:
            raise NotImplementedError(f'spin type `{spin_type}` is not supported.')

        if initial_magnetic_moments is not None and spin_type is not SpinType.COLLINEAR:
            raise ValueError(f'`initial_magnetic_moments` is specified but spin type `{spin_type}` is incompatible.')

        inputs = cls.get_protocol_inputs(protocol, overrides)

        meta_parameters = inputs.pop('meta_parameters')
        pseudo_family = inputs.pop('pseudo_family')

        natoms = len(structure.sites)

        try:
            pseudo_set = (PseudoDojoFamily, SsspFamily, CutoffsPseudoPotentialFamily)
            pseudo_family = orm.QueryBuilder().append(pseudo_set, filters={'label': pseudo_family}).one()[0]
        except exceptions.NotExistent as exception:
            raise ValueError(
                f'required pseudo family `{pseudo_family}` is not installed. Please use `aiida-pseudo install` to'
                'install it.'
            ) from exception

        try:
            cutoff_wfc, cutoff_rho = pseudo_family.get_recommended_cutoffs(structure=structure, unit='Ry')
            pseudos = pseudo_family.get_pseudos(structure=structure)
        except ValueError as exception:
            raise ValueError(
                f'failed to obtain recommended cutoffs for pseudo family `{pseudo_family}`: {exception}'
            ) from exception

        # Update the parameters based on the protocol inputs
        parameters = inputs['pw']['parameters']
        parameters['CONTROL']['etot_conv_thr'] = natoms * meta_parameters['etot_conv_thr_per_atom']
        parameters['ELECTRONS']['conv_thr'] = natoms * meta_parameters['conv_thr_per_atom']
        parameters['SYSTEM']['ecutwfc'] = cutoff_wfc
        parameters['SYSTEM']['ecutrho'] = cutoff_rho

        if electronic_type is ElectronicType.INSULATOR:
            parameters['SYSTEM']['occupations'] = 'fixed'
            parameters['SYSTEM'].pop('degauss')
            parameters['SYSTEM'].pop('smearing')

        if spin_type is SpinType.COLLINEAR:
            starting_magnetization = get_starting_magnetization(structure, pseudo_family, initial_magnetic_moments)
            parameters['SYSTEM']['starting_magnetization'] = starting_magnetization
            parameters['SYSTEM']['nspin'] = 2

        # If overrides are provided, they are considered absolute
        if overrides:
            parameter_overrides = overrides.get('pw', {}).get('parameters', {})
            parameters = recursive_merge(parameters, parameter_overrides)

            pseudos_overrides = overrides.get('pw', {}).get('pseudos', {})
            pseudos = recursive_merge(pseudos, pseudos_overrides)

        metadata = inputs['pw']['metadata']

        if options:
            metadata['options'] = recursive_merge(inputs['pw']['metadata']['options'], options)

        # pylint: disable=no-member
        builder = cls.get_builder()
        builder.pw['code'] = code
        builder.pw['pseudos'] = pseudos
        builder.pw['structure'] = structure
        builder.pw['parameters'] = orm.Dict(parameters)
        builder.pw['metadata'] = metadata
        if 'settings' in inputs['pw']:
            builder.pw['settings'] = orm.Dict(inputs['pw']['settings'])
        if 'parallelization' in inputs['pw']:
            builder.pw['parallelization'] = orm.Dict(inputs['pw']['parallelization'])
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])
        if 'kpoints' in inputs:
            builder.kpoints = inputs['kpoints']
        else:
            builder.kpoints_distance = orm.Float(inputs['kpoints_distance'])
        builder.kpoints_force_parity = orm.Bool(inputs['kpoints_force_parity'])
        # pylint: enable=no-member

        return builder

    def setup(self):
        """Call the ``setup`` of the ``BaseRestartWorkChain`` and create the inputs dictionary in ``self.ctx.inputs``.

        This ``self.ctx.inputs`` dictionary will be used by the ``BaseRestartWorkChain`` to submit the calculations
        in the internal loop.

        The ``parameters`` and ``settings`` input ``Dict`` nodes are converted into a regular dictionary and the
        default namelists for the ``parameters`` are set to empty dictionaries if not specified.
        """
        super().setup()
        self.ctx.inputs = AttributeDict(self.exposed_inputs(PwCalculation, 'pw'))

        self.ctx.inputs.parameters = self.ctx.inputs.parameters.get_dict()
        self.ctx.inputs.parameters.setdefault('CONTROL', {})
        self.ctx.inputs.parameters.setdefault('ELECTRONS', {})
        self.ctx.inputs.parameters.setdefault('SYSTEM', {})
        self.ctx.inputs.parameters['CONTROL'].setdefault('calculation', 'scf')

        self.ctx.inputs.settings = self.ctx.inputs.settings.get_dict() if 'settings' in self.ctx.inputs else {}

        # If a ``parent_folder`` is specified, automatically set the parameters for a ``RestartType.Full`` unless the
        # ``CONTROL.restart_mode`` has explicitly been set to ``from_scratch``. In that case, the user most likely set
        # that, and we do not want to override it.
        restart_mode = self.ctx.inputs.parameters['CONTROL'].get('restart_mode', None)
        if 'parent_folder' in self.ctx.inputs and restart_mode != 'from_scratch':
            self.set_restart_type(RestartType.FULL, self.ctx.inputs.parent_folder)

    def validate_kpoints(self):
        """Validate the inputs related to k-points.

        Either an explicit `KpointsData` with given mesh/path, or a desired k-points distance should be specified. In
        the case of the latter, the `KpointsData` will be constructed for the input `StructureData` using the
        `create_kpoints_from_distance` calculation function.
        """
        if all(key not in self.inputs for key in ['kpoints', 'kpoints_distance']):
            return self.exit_codes.ERROR_INVALID_INPUT_KPOINTS

        try:
            kpoints = self.inputs.kpoints
        except AttributeError:
            inputs = {
                'structure': self.inputs.pw.structure,
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
        max_seconds_factor = self.defaults.delta_factor_max_seconds
        max_seconds = max_wallclock_seconds * max_seconds_factor
        self.ctx.inputs.parameters['CONTROL']['max_seconds'] = max_seconds

    def set_restart_type(self, restart_type, parent_folder=None):
        """Set the restart type for the next iteration."""

        if parent_folder is None and restart_type != RestartType.FROM_SCRATCH:
            raise ValueError('When not restarting from scratch, a `parent_folder` must be provided.')

        if restart_type == RestartType.FROM_SCRATCH:
            self.ctx.inputs.parameters['CONTROL']['restart_mode'] = 'from_scratch'
            self.ctx.inputs.parameters['ELECTRONS'].pop('startingpot', None)
            self.ctx.inputs.parameters['ELECTRONS'].pop('startingwfc', None)
            self.ctx.inputs.pop('parent_folder', None)

        elif restart_type == RestartType.FULL:
            self.ctx.inputs.parameters['CONTROL']['restart_mode'] = 'restart'
            self.ctx.inputs.parameters['ELECTRONS'].pop('startingpot', None)
            self.ctx.inputs.parameters['ELECTRONS'].pop('startingwfc', None)
            self.ctx.inputs.parent_folder = parent_folder

        elif restart_type == RestartType.FROM_CHARGE_DENSITY:
            self.ctx.inputs.parameters['CONTROL']['restart_mode'] = 'from_scratch'
            self.ctx.inputs.parameters['ELECTRONS']['startingpot'] = 'file'
            self.ctx.inputs.parameters['ELECTRONS'].pop('startingwfc', None)
            self.ctx.inputs.parent_folder = parent_folder

        elif restart_type == RestartType.FROM_WAVE_FUNCTIONS:
            self.ctx.inputs.parameters['CONTROL']['restart_mode'] = 'from_scratch'
            self.ctx.inputs.parameters['ELECTRONS'].pop('startingpot', None)
            self.ctx.inputs.parameters['ELECTRONS']['startingwfc'] = 'file'
            self.ctx.inputs.parent_folder = parent_folder

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
            self.report(f'required key "{key}" in automatic_parallelization input not found')
            return self.exit_codes.ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_MISSING_KEY

        if remaining_keys:
            self.report(
                f"detected unrecognized keys in the automatic_parallelization input: {' '.join(remaining_keys)}"
            )
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

        self.report(f'launching initialization {running.pk}<{self.ctx.process_name}>')

        return ToContext(calculation_init=running)

    def inspect_init(self):
        """Use the initialization `PwCalculation` to determine the required resource and parallelization settings."""
        calculation = self.ctx.calculation_init

        if not calculation.is_finished_ok:
            return self.exit_codes.ERROR_INITIALIZATION_CALCULATION_FAILED

        # Get automated parallelization settings
        parallelization = get_pw_parallelization_parameters(calculation, **self.ctx.automatic_parallelization)

        # Note: don't do this at home, we are losing provenance here. This should be done by a calculation function
        node = orm.Dict(parallelization).store()
        self.out('automatic_parallelization', node)
        self.report(f'results of automatic parallelization in {node.__class__.__name__}<{node.pk}>')

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

    def prepare_process(self):
        """Prepare the inputs for the next calculation."""
        max_wallclock_seconds = self.ctx.inputs.metadata.options.get('max_wallclock_seconds', None)

        if max_wallclock_seconds is not None and 'max_seconds' not in self.ctx.inputs.parameters['CONTROL']:
            self.set_max_seconds(max_wallclock_seconds)

    def report_error_handled(self, calculation, action):
        """Report an action taken for a calculation that has failed.

        This should be called in a registered error handler if its condition is met and an action was taken.

        :param calculation: the failed calculation node
        :param action: a string message with the action taken
        """
        arguments = [calculation.process_label, calculation.pk, calculation.exit_status, calculation.exit_message]
        self.report('{}<{}> failed with exit status {}: {}'.format(*arguments))
        self.report(f'Action taken: {action}')

    @process_handler(exit_codes=ExitCode(0))
    def sanity_check_insufficient_bands(self, calculation):
        """Perform a sanity check on the band occupations of a  successfully converged calculation.

        Verify that the occupation of the last band is below a certain threshold, unless `occupations` was explicitly
        set to `fixed` in the input parameters. If this is violated, the calculation used too few bands and cannot be
        trusted. The number of bands is increased and the calculation is restarted, using the charge density from the
        previous calculation.
        """
        from aiida_quantumespresso.utils.bands import get_highest_occupied_band

        occupations = calculation.inputs.parameters.base.attributes.get('SYSTEM', {}).get('occupations', None)

        if occupations is None:
            self.report(
                '`SYSTEM.occupations` parameter is not defined: performing band occupation check. '
                'If you want to disable this, explicitly set `SYSTEM.occupations` to `fixed`.'
            )

        # Only skip the check on the highest band occupation if `occupations` was explicitly set to `fixed`.
        if occupations == 'fixed':
            return

        try:
            bands = calculation.outputs.output_band
        except AttributeError:
            args = [self.ctx.process_name, calculation.pk]
            self.report('{}<{}> does not have `output_band` output, skipping sanity check.'.format(*args))
            return

        try:
            get_highest_occupied_band(bands)
        except ValueError as exception:
            args = [self.ctx.process_name, calculation.pk]
            self.report('{}<{}> run with smearing and highest band is occupied'.format(*args))
            self.report(f'BandsData<{bands.pk}> has invalid occupations: {exception}')
            self.report(f'{calculation.process_label}<{calculation.pk}> had insufficient bands')

            nbnd_cur = calculation.outputs.output_parameters.get_dict()['number_of_bands']
            nbnd_new = nbnd_cur + max(int(nbnd_cur * self.defaults.delta_factor_nbnd), self.defaults.delta_minimum_nbnd)
            self.ctx.inputs.parameters['SYSTEM']['nbnd'] = nbnd_new

            self.set_restart_type(RestartType.FROM_CHARGE_DENSITY, calculation.outputs.remote_folder)
            self.report(
                f'Action taken: increased number of bands to {nbnd_new} and restarting from the previous charge '
                'density.'
            )

            return ProcessHandlerReport(True)

    @process_handler(priority=600)
    def handle_unrecoverable_failure(self, calculation):
        """Handle calculations with an exit status below 400 which are unrecoverable, so abort the work chain."""
        if calculation.is_failed and calculation.exit_status < 400:
            self.report_error_handled(calculation, 'unrecoverable error, aborting...')
            return ProcessHandlerReport(True, self.exit_codes.ERROR_UNRECOVERABLE_FAILURE)

    @process_handler(priority=590, exit_codes=[])
    def handle_known_unrecoverable_failure(self, calculation):
        """Handle calculations with an exit status that correspond to a known failure mode that are unrecoverable.

        These failures may always be unrecoverable or at some point a handler may be devised.
        """
        self.report_error_handled(calculation, 'known unrecoverable failure detected, aborting...')
        return ProcessHandlerReport(True, self.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE)

    @process_handler(
        priority=585,
        exit_codes=[
            PwCalculation.exit_codes.ERROR_COMPUTING_CHOLESKY,
            PwCalculation.exit_codes.ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED,
        ]
    )
    def handle_diagonalization_errors(self, calculation):
        """Handle known issues related to the diagonalization.

        Switch to ``diagonalization = 'cg'`` if not already running with this setting, and restart from the charge
        density. In case the run already used conjugate gradient diagonalization, abort.
        """
        if self.ctx.inputs.parameters['ELECTRONS'].get('diagonalization', None) == 'cg':
            action = 'found diagonalization issues but already running with conjugate gradient algorithm, aborting...'
            self.report_error_handled(calculation, action)
            return ProcessHandlerReport(True, self.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE)

        self.ctx.inputs.parameters['ELECTRONS']['diagonalization'] = 'cg'
        action = 'found diagonalization issues, switching to conjugate gradient diagonalization.'
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True)

    @process_handler(priority=580, exit_codes=[
        PwCalculation.exit_codes.ERROR_OUT_OF_WALLTIME,
    ])
    def handle_out_of_walltime(self, calculation):
        """Handle `ERROR_OUT_OF_WALLTIME` exit code.

        In this case the calculation shut down neatly and we can simply restart. We consider two cases:

        1. If the structure is unchanged, we do a full restart.
        2. If the structure has changed during the calculation, we restart from scratch.
        """
        try:
            self.ctx.inputs.structure = calculation.outputs.output_structure
        except exceptions.NotExistent:
            self.set_restart_type(RestartType.FULL, calculation.outputs.remote_folder)
            self.report_error_handled(calculation, 'simply restart from the last calculation')
        else:
            self.set_restart_type(RestartType.FROM_SCRATCH)
            self.report_error_handled(calculation, 'out of walltime: structure changed so restarting from scratch')

        return ProcessHandlerReport(True)

    @process_handler(
        priority=570, exit_codes=[
            PwCalculation.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF,
        ]
    )
    def handle_vcrelax_converged_except_final_scf(self, calculation):
        """Handle `ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF` exit code.

        Convergence reached in `vc-relax` except thresholds exceeded in final scf: consider as converged.
        """
        self.ctx.is_finished = True
        action = 'ionic convergence thresholds met except in final scf: consider structure relaxed.'
        self.report_error_handled(calculation, action)
        self.results()  # Call the results method to attach the output nodes
        return ProcessHandlerReport(True, self.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF)

    @process_handler(
        priority=560,
        exit_codes=[
            PwCalculation.exit_codes.ERROR_IONIC_CONVERGENCE_NOT_REACHED,
            PwCalculation.exit_codes.ERROR_IONIC_CYCLE_EXCEEDED_NSTEP,
            PwCalculation.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE,
            PwCalculation.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE,
        ]
    )
    def handle_relax_recoverable_ionic_convergence_error(self, calculation):
        """Handle various exit codes for recoverable `relax` calculations with failed ionic convergence.

        These exit codes signify that the ionic convergence thresholds were not met, but the output structure is usable,
        so the solution is to simply restart from scratch but from the output structure.
        """
        self.ctx.inputs.structure = calculation.outputs.output_structure
        action = 'no ionic convergence but clean shutdown: restarting from scratch but using output structure.'

        self.set_restart_type(RestartType.FROM_SCRATCH)
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True)

    @process_handler(
        priority=550,
        exit_codes=[
            PwCalculation.exit_codes.ERROR_IONIC_CYCLE_ELECTRONIC_CONVERGENCE_NOT_REACHED,
            PwCalculation.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_FINAL_SCF_FAILED,
        ]
    )
    def handle_relax_recoverable_electronic_convergence_error(self, calculation):
        """Handle various exit codes for recoverable `relax` calculations with failed electronic convergence.

        These exit codes signify that the electronic convergence thresholds were not met, but the output structure is
        usable, so the solution is to simply restart from scratch but from the output structure and with a reduced
        ``mixing_beta``.
        """
        factor = self.defaults.delta_factor_mixing_beta
        mixing_beta = self.ctx.inputs.parameters.get('ELECTRONS', {}).get('mixing_beta', self.defaults.qe.mixing_beta)
        mixing_beta_new = mixing_beta * factor

        self.ctx.inputs.parameters['ELECTRONS']['mixing_beta'] = mixing_beta_new
        self.ctx.inputs.structure = calculation.outputs.output_structure
        action = (
            f'no electronic convergence but clean shutdown: reduced beta mixing from {mixing_beta} to {mixing_beta_new}'
            'restarting from scratch but using output structure.'
        )

        self.set_restart_type(RestartType.FROM_SCRATCH)
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True)

    @process_handler(priority=410, exit_codes=[
        PwCalculation.exit_codes.ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED,
    ])
    def handle_electronic_convergence_not_reached(self, calculation):
        """Handle `ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED` error.

        Decrease the mixing beta and fully restart from the previous calculation.
        """
        factor = self.defaults.delta_factor_mixing_beta
        mixing_beta = self.ctx.inputs.parameters.get('ELECTRONS', {}).get('mixing_beta', self.defaults.qe.mixing_beta)
        mixing_beta_new = mixing_beta * factor

        self.ctx.inputs.parameters['ELECTRONS']['mixing_beta'] = mixing_beta_new
        action = f'reduced beta mixing from {mixing_beta} to {mixing_beta_new} and restarting from the last calculation'

        self.set_restart_type(RestartType.FULL, calculation.outputs.remote_folder)
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True)

    @process_handler(priority=420, exit_codes=[
        PwCalculation.exit_codes.WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED,
    ])
    def handle_electronic_convergence_warning(self, calculation):
        """Handle `WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED`: consider finished."""
        self.ctx.is_finished = True
        action = 'electronic convergence not reached but inputs say this is ok: consider finished.'
        self.report_error_handled(calculation, action)
        self.results()  # Call the results method to attach the output nodes
        return ProcessHandlerReport(True, self.exit_codes.WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED)
