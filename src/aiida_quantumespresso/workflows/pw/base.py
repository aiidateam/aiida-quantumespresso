"""Workchain to run a Quantum ESPRESSO pw.x calculation with automated error handling and restarts."""

import warnings

from aiida import orm
from aiida.common import AttributeDict, exceptions
from aiida.common.lang import type_check
from aiida.engine import (
    BaseRestartWorkChain,
    ExitCode,
    ProcessHandlerReport,
    process_handler,
    while_,
)
from aiida.plugins import CalculationFactory, GroupFactory
from aiida_pseudo.groups.mixins import RecommendedCutoffMixin

from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import (
    create_kpoints_from_distance,
)
from aiida_quantumespresso.calculations.functions.rattle_structure import rattle_structure
from aiida_quantumespresso.common.types import ElectronicType, RestartType, SpinType
from aiida_quantumespresso.tools.monitors.accuracy_stuck import ACCURACY_STUCK_MESSAGE_PREFIX
from aiida_quantumespresso.utils.defaults.calculation import pw as qe_defaults

from ..protocols.utils import ProtocolMixin

PwCalculation = CalculationFactory('quantumespresso.pw')
PseudoPotentialFamily = GroupFactory('pseudo.family')


class PwBaseWorkChain(ProtocolMixin, BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO pw.x calculation with automated error handling and restarts."""

    _process_class = PwCalculation

    defaults = AttributeDict(
        {
            'qe': qe_defaults,
            'delta_threshold_degauss': 30,
            'delta_factor_degauss': 0.1,
            'delta_factor_mixing_beta': 0.8,
            'delta_factor_max_seconds': 0.95,
            'delta_factor_nbnd': 0.05,
            'delta_minimum_nbnd': 4,
            'delta_factor_trust_radius_min': 0.1,
            'delta_factor_upscale': 0.3,
            'rattle_stdev': 0.01,
            'rattle_symmetry_threshold': 5,
        }
    )

    @classmethod
    def define(cls, spec):
        """Define the process specification."""

        super().define(spec)
        spec.expose_inputs(PwCalculation, namespace='pw', exclude=('kpoints',))
        spec.input(
            'kpoints',
            valid_type=orm.KpointsData,
            required=False,
            help='An explicit k-points list or mesh. Either this or `kpoints_distance` has to be provided.',
        )
        spec.input(
            'kpoints_distance',
            valid_type=orm.Float,
            required=False,
            help='The minimum desired distance in 1/Å between k-points in reciprocal space. The explicit k-points will '
            'be generated automatically by a calculation function based on the input structure.',
        )
        spec.input(
            'kpoints_force_parity',
            valid_type=orm.Bool,
            required=False,
            help='Optional input when constructing the k-points based on a desired `kpoints_distance`. Setting this to '
            '`True` will force the k-point mesh to have an even number of points along each lattice vector except '
            'for any non-periodic directions.',
        )

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

        spec.expose_outputs(PwCalculation)

        spec.exit_code(
            201,
            'ERROR_INVALID_INPUT_PSEUDO_POTENTIALS',
            message='The explicit `pseudos` or `pseudo_family` could not be used to get the necessary pseudos.',
        )
        spec.exit_code(
            202,
            'ERROR_INVALID_INPUT_KPOINTS',
            message='Neither the `kpoints` nor the `kpoints_distance` input was specified.',
        )
        spec.exit_code(
            203,
            'ERROR_INVALID_INPUT_RESOURCES',
            message='Neither the `options` nor `automatic_parallelization` input was specified. '
            'This exit status has been deprecated as the check it corresponded to was incorrect.',
        )
        spec.exit_code(
            204,
            'ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED',
            message='The `metadata.options` did not specify both `resources.num_machines` and `max_wallclock_seconds`. '
            'This exit status has been deprecated as the check it corresponded to was incorrect.',
        )
        spec.exit_code(
            210,
            'ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_MISSING_KEY',
            message='Required key for `automatic_parallelization` was not specified.'
            'This exit status has been deprecated as the automatic parallellization feature was removed.',
        )
        spec.exit_code(
            211,
            'ERROR_INVALID_INPUT_AUTOMATIC_PARALLELIZATION_UNRECOGNIZED_KEY',
            message='Unrecognized keys were specified for `automatic_parallelization`.'
            'This exit status has been deprecated as the automatic parallellization feature was removed.',
        )
        spec.exit_code(
            300,
            'ERROR_UNRECOVERABLE_FAILURE',
            message='[deprecated] The calculation failed with an unidentified unrecoverable error.',
        )
        spec.exit_code(
            310, 'ERROR_KNOWN_UNRECOVERABLE_FAILURE', message='The calculation failed with a known unrecoverable error.'
        )
        spec.exit_code(320, 'ERROR_INITIALIZATION_CALCULATION_FAILED', message='The initialization calculation failed.')
        spec.exit_code(
            501,
            'ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF',
            message='Then ionic minimization cycle converged but the thresholds are exceeded in the final SCF.',
        )
        spec.exit_code(
            710,
            'WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED',
            message='The electronic minimization cycle did not reach self-consistency, but `scf_must_converge` '
            'is `False` and/or `electron_maxstep` is 0.',
        )

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
        **_,
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
        from aiida_quantumespresso.workflows.protocols.utils import (
            get_magnetization,
            recursive_merge,
        )

        if isinstance(code, str):
            code = orm.load_code(code)

        type_check(code, orm.AbstractCode)
        type_check(electronic_type, ElectronicType)
        type_check(spin_type, SpinType)

        if electronic_type not in [ElectronicType.METAL, ElectronicType.INSULATOR]:
            raise NotImplementedError(f'electronic type `{electronic_type}` is not supported.')

        if initial_magnetic_moments is not None and spin_type == SpinType.NONE:
            raise ValueError(f'`initial_magnetic_moments` is specified but spin type `{spin_type}` is incompatible.')

        inputs = cls.get_protocol_inputs(protocol, overrides)

        meta_parameters = inputs.pop('meta_parameters')
        pseudo_family = inputs.pop('pseudo_family')

        if spin_type is SpinType.SPIN_ORBIT and 'pseudo_family' not in (overrides or {}):
            pseudo_family = 'PseudoDojo/0.4/PBEsol/FR/standard/upf'

        natoms = len(structure.sites)

        # Update the parameters based on the protocol inputs
        parameters = inputs['pw']['parameters']

        system_overrides = (overrides or {}).get('pw', {}).get('parameters', {}).get('SYSTEM', {})
        cutoffs_in_overrides = all(key in system_overrides for key in ('ecutwfc', 'ecutrho'))

        if overrides and 'pseudos' in overrides.get('pw', {}):
            pseudos = overrides['pw']['pseudos']

            if sorted(pseudos.keys()) != sorted(structure.get_kind_names()):
                raise ValueError(f'`pseudos` override needs one value for each of the {len(structure.kinds)} kinds.')

            if not cutoffs_in_overrides:
                raise ValueError(
                    'When overriding the pseudo potentials, both `ecutwfc` and `ecutrho` cutoffs should be '
                    f'provided in the `overrides`: {overrides}'
                )

        else:
            query = orm.QueryBuilder().append(PseudoPotentialFamily, filters={'label': pseudo_family})

            try:
                pseudo_family = query.one()[0]
            except exceptions.NotExistent as exception:
                raise ValueError(
                    f'required pseudo family `{pseudo_family}` is not installed. Please use `aiida-pseudo install` to'
                    'install it.'
                ) from exception
            except exceptions.MultipleObjectsError as exception:
                matches = ', '.join(f'`{family}`' for family in query.all(flat=True))
                raise ValueError(
                    f'the label `{pseudo_family}` matches more than one installed pseudo family: {matches}. Please '
                    'delete or relabel all but one, or pass the `pseudos` in the `overrides` instead.'
                ) from exception

            # Families that do not define recommended cutoffs can still be used, as long as the `overrides` provide
            # both cutoffs themselves, since these are applied further down and take precedence anyway.
            if not cutoffs_in_overrides:
                if not isinstance(pseudo_family, RecommendedCutoffMixin):
                    raise ValueError(
                        f'pseudo family `{pseudo_family}` cannot recommend cutoffs. Provide both `ecutwfc` and '
                        '`ecutrho` in the `overrides`, or use a family that recommends them.'
                    )

                try:
                    parameters['SYSTEM']['ecutwfc'], parameters['SYSTEM']['ecutrho'] = (
                        pseudo_family.get_recommended_cutoffs(structure=structure, unit='Ry')
                    )
                except ValueError as exception:
                    raise ValueError(
                        f'failed to obtain recommended cutoffs for pseudo family `{pseudo_family}`: {exception} '
                        'If the family does not define any, specify both `ecutwfc` and `ecutrho` in the `overrides`.'
                    ) from exception

            pseudos = pseudo_family.get_pseudos(structure=structure)

        parameters['CONTROL']['etot_conv_thr'] = natoms * meta_parameters['etot_conv_thr_per_atom']
        parameters['ELECTRONS']['conv_thr'] = natoms * meta_parameters['conv_thr_per_atom']

        pbc_parameter_overrides = {
            (True, True, False): {'SYSTEM': {'assume_isolated': '2D'}},
        }
        if structure.pbc != (True, True, True):
            # 0D, 1D, or 2D
            if structure.pbc.count(True) == 2 and structure.pbc not in pbc_parameter_overrides:
                raise ValueError(f'2D-periodic structures must be periodic in the x-y plane, got `{structure.pbc}`.')
            warnings.warn(
                'This protocol was developed for fully periodic (i.e. 3D) systems. Use `overrides` to provide '
                'any relevant keywords for handling aperiodicity, and proceed with caution.'
            )
        parameters = recursive_merge(parameters, pbc_parameter_overrides.get(structure.pbc, {}))

        if electronic_type is ElectronicType.INSULATOR:
            parameters['SYSTEM']['occupations'] = 'fixed'
            parameters['SYSTEM'].pop('degauss')
            parameters['SYSTEM'].pop('smearing')

        if spin_type in [SpinType.COLLINEAR, SpinType.SPIN_ORBIT, SpinType.NON_COLLINEAR]:
            magnetization = get_magnetization(
                structure=structure,
                z_valences={kind.name: pseudos[kind.name].z_valence for kind in structure.kinds},
                initial_magnetic_moments=initial_magnetic_moments,
                spin_type=spin_type,
            )
            parameters['SYSTEM']['starting_magnetization'] = magnetization['starting_magnetization']

            if spin_type is SpinType.COLLINEAR:
                parameters['SYSTEM']['nspin'] = 2

            if spin_type in [SpinType.SPIN_ORBIT, SpinType.NON_COLLINEAR]:
                parameters['SYSTEM']['angle1'] = magnetization['angle1']
                parameters['SYSTEM']['angle2'] = magnetization['angle2']
                parameters['SYSTEM']['noncolin'] = True
                parameters['SYSTEM']['nspin'] = 4

            if spin_type is SpinType.SPIN_ORBIT:
                parameters['SYSTEM']['lspinorb'] = True

        # If overrides are provided, they are considered absolute
        if overrides:
            parameter_overrides = overrides.get('pw', {}).get('parameters', {})
            parameters = recursive_merge(parameters, parameter_overrides)

            # if tot_magnetization in overrides , remove starting_magnetization from parameters
            if parameters.get('SYSTEM', {}).get('tot_magnetization') is not None:
                parameters.setdefault('SYSTEM', {}).pop('starting_magnetization', None)

        metadata = inputs['pw']['metadata']

        if options:
            metadata['options'] = recursive_merge(metadata['options'], options)

        metadata['options'] = cls.set_default_resources(metadata['options'], code.computer.scheduler_type)

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
        builder.max_iterations = orm.Int(inputs['max_iterations'])

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

        calculation_type = self.ctx.inputs.parameters['CONTROL'].get('calculation', None)
        if calculation_type in ['relax', 'md']:
            self.ctx.inputs.parameters.setdefault('IONS', {})
        if calculation_type in ['vc-relax', 'vc-md']:
            self.ctx.inputs.parameters.setdefault('IONS', {})
            self.ctx.inputs.parameters.setdefault('CELL', {})

        self.ctx.inputs.settings = self.ctx.inputs.settings.get_dict() if 'settings' in self.ctx.inputs else {}

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
                'metadata': {'call_link_label': 'create_kpoints_from_distance'},
            }
            kpoints = create_kpoints_from_distance(**inputs)

        self.ctx.inputs.kpoints = kpoints

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

    def prepare_process(self):
        """Prepare the inputs for the next calculation."""
        max_wallclock_seconds = self.ctx.inputs.metadata.options.get('max_wallclock_seconds', None)

        if max_wallclock_seconds is not None and 'max_seconds' not in self.ctx.inputs.parameters['CONTROL']:
            max_seconds = max_wallclock_seconds * self.defaults.delta_factor_max_seconds
            self.ctx.inputs.parameters['CONTROL']['max_seconds'] = max_seconds

    def report_error_handled(self, calculation, action):
        """Report an action taken for a calculation that has failed.

        This should be called in a registered error handler if its condition is met and an action was taken.

        :param calculation: the failed calculation node
        :param action: a string message with the action taken
        """
        arguments = [
            calculation.process_label,
            calculation.pk,
            calculation.exit_status,
            calculation.exit_message,
        ]
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
            nbnd_new = nbnd_cur + max(
                int(nbnd_cur * self.defaults.delta_factor_nbnd),
                self.defaults.delta_minimum_nbnd,
            )
            self.ctx.inputs.parameters['SYSTEM']['nbnd'] = nbnd_new

            self.set_restart_type(RestartType.FROM_CHARGE_DENSITY, calculation.outputs.remote_folder)
            self.report(
                f'Action taken: increased number of bands to {nbnd_new} and restarting from the previous charge '
                'density.'
            )

            return ProcessHandlerReport(True)

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
            PwCalculation.exit_codes.ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE,
            PwCalculation.exit_codes.ERROR_ZHEGVD_FAILED,
            PwCalculation.exit_codes.ERROR_QR_FAILED,
            PwCalculation.exit_codes.ERROR_EIGENVECTOR_CONVERGENCE,
            PwCalculation.exit_codes.ERROR_BROYDEN_FACTORIZATION,
        ],
    )
    def handle_diagonalization_errors(self, calculation):
        """Handle known issues related to the diagonalization.

        We use the following strategy. When a diagonalization algorithm fails, we try using an other one
        still not used. Conjugate gradient (CG) is kept as last option, as it is the slowest among the
        available ones, but on the contrary it is the most stable as well, thus kept as `last resort`.

        Once the error handler has tried all ``diagonalization`` options, abort.
        """
        current = self.ctx.inputs.parameters['ELECTRONS'].get('diagonalization', 'david')

        if 'diagonalizations' not in self.ctx:
            # Initialize a list to track diagonalisations that haven't been tried in reverse order or preference
            self.ctx.diagonalizations = [value for value in ['cg', 'paro', 'david'] if value != current.lower()]

        try:
            new = self.ctx.diagonalizations.pop()
            self.ctx.inputs.parameters['ELECTRONS']['diagonalization'] = new
            action = f'found diagonalization issues for ``{current}``, switching to ``{new}`` diagonalization.'
            self.report_error_handled(calculation, action)
            return ProcessHandlerReport(True)
        except IndexError:
            action = 'found diagonalization issues but already exploited all supported algorithms, aborting...'
            self.report_error_handled(calculation, action)
            return ProcessHandlerReport(True, self.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE)

    @process_handler(
        priority=580,
        exit_codes=[
            PwCalculation.exit_codes.ERROR_OUT_OF_WALLTIME,
        ],
    )
    def handle_out_of_walltime(self, calculation):
        """Handle `ERROR_OUT_OF_WALLTIME` exit code.

        In this case the calculation shut down cleanly and we can do a full restart.
        """
        if 'output_structure' in calculation.outputs:
            self.ctx.inputs.structure = calculation.outputs.output_structure

        self.set_restart_type(RestartType.FULL, calculation.outputs.remote_folder)
        self.report_error_handled(calculation, "restarting in full with `CONTROL.restart_mode` = 'restart'")

        return ProcessHandlerReport(True)

    @process_handler(
        priority=575,
        exit_codes=[
            PwCalculation.exit_codes.ERROR_IONIC_INTERRUPTED_PARTIAL_TRAJECTORY,
        ],
    )
    def handle_ionic_interrupted_partial_trajectory(self, calculation):
        """Handle `ERROR_IONIC_INTERRUPTED_PARTIAL_TRAJECTORY` exit code.

        In this case the calculation got interrupted during an ionic optimization due to a problem that is likely
        transient, so we can restart from the last output structure. Note that since the job got interrupted the charge
        density and wave functions are likely corrupt so those cannot be used in the restart.
        """
        self.ctx.inputs.structure = calculation.outputs.output_structure
        self.set_restart_type(RestartType.FROM_SCRATCH)
        self.report_error_handled(calculation, 'restarting from scratch from the last output structure')
        return ProcessHandlerReport(True)

    @process_handler(
        priority=570,
        exit_codes=[
            PwCalculation.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF,
        ],
    )
    def handle_vcrelax_converged_except_final_scf(self, calculation):
        """Handle `ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF` exit code.

        Convergence reached in `vc-relax` except thresholds exceeded in final scf: consider as converged.
        """
        self.ctx.is_finished = True
        action = 'consider structure final, but report exit code.'
        self.report_error_handled(calculation, action)
        self.results()  # Call the results method to attach the output nodes
        return ProcessHandlerReport(True, self.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF)

    @process_handler(
        priority=602,
        exit_codes=[
            PwCalculation.exit_codes.STOPPED_BY_MONITOR,
        ],
    )
    def handle_electronic_convergence_stuck(self, calculation):
        """Handle `STOPPED_BY_MONITOR` error.

        This is needed when the scf convergence is stuck during a relax calculation, i.e. when the algorithm in
        Quantum ESPRESSO tries to decrease the `conv_thr` during the ionic minimization, to ensure that forces are
        small enough. To solve it, we restart by reducing the `IONS.upscale` parameter (by default, it is 100).
        """
        if calculation.exit_message is None or not calculation.exit_message.startswith(ACCURACY_STUCK_MESSAGE_PREFIX):
            # The exit code is shared by all monitors, so the message identifies the accuracy-stuck monitor. Any
            # other message comes from an unrelated monitor: return `None` so the failure counts as unhandled.
            return None

        calculation_type = self.ctx.inputs.parameters['CONTROL'].get('calculation', 'scf')

        if calculation_type in ('relax', 'vc-relax'):
            ions_namelist = self.ctx.inputs.parameters.get('IONS', {})
            upscale = ions_namelist.get('upscale', qe_defaults.upscale)

            if upscale <= 1:
                self.report_error_handled(
                    calculation, 'electronic convergence stuck but `upscale` already at its minimum, aborting...'
                )
                return ProcessHandlerReport(True, self.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE)

            new_upscale = max(int(upscale * self.defaults.delta_factor_upscale), 1)
            ions_namelist['upscale'] = new_upscale
            self.ctx.inputs.parameters['IONS'] = ions_namelist
            action = f'reduced upscale from {upscale} to {new_upscale} and restarting from the last calculation'
        else:
            action = 'restarting from the last calculation'

        self.set_restart_type(RestartType.FULL, calculation.outputs.remote_folder)
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True)

    @process_handler(
        priority=601,
        exit_codes=[
            PwCalculation.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE,
            PwCalculation.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE,
        ],
    )
    def handle_relax_recoverable_ionic_convergence_bfgs_history_error(self, calculation):
        """Handle failure of the ionic minimization algorithm (BFGS).

        When BFGS history fails, this can mean two things: the structure is close to the global minimum, but the
        moves the algorithm wants to do are smaller than `trust_radius_min`, or the structure is close to a local
        minimum (hard to detect). For the first case, one can lower the trust radius; for the second one, one can
        exploit a different algorithm, e.g. `damp` (and `damp-w` for vc-relax). Additionally, in the case of a `relax`
        calculation with few symmetries, the structure is rattled once (with a 0.01 Angstrom standard deviation) to
        escape from a possible (hard) local minimum. A second occurrence of the same failure after the rattle is
        considered unrecoverable.
        """
        if not hasattr(self.ctx, 'rattled'):
            self.ctx.rattled = False

        trust_radius_min = self.ctx.inputs.parameters['IONS'].get('trust_radius_min', qe_defaults.trust_radius_min)
        calculation_type = self.ctx.inputs.parameters['CONTROL'].get('calculation', 'relax')

        if calculation_type == 'relax':
            n_symmetries = calculation.outputs.output_parameters.get_dict().get('number_of_symmetries')
            if n_symmetries is not None and n_symmetries < self.defaults.rattle_symmetry_threshold:
                # Handle the case of a 'hard' local minimum by rattling the structure once, then a second
                # occurrence is considered unrecoverable. Restore the original upscale value (i.e. from the
                # initial input parameters, not from the possibly reduced value set by the handler dealing with
                # the `electronic_convergence_stuck` error), since the structure is no longer close to the
                # global minimum. Otherwise, the reduced upscale would make it impossible to reach convergence.
                if not self.ctx.rattled:
                    self.ctx.rattled = True
                    self.ctx.inputs.structure = rattle_structure(
                        calculation.outputs.output_structure, orm.Float(self.defaults.rattle_stdev)
                    )
                    self.ctx.inputs.parameters['IONS']['upscale'] = (
                        self.inputs.pw.parameters.get_dict().get('IONS', {}).get('upscale', qe_defaults.upscale)
                    )
                    action = (
                        'bfgs history (ionic only) failure: rattling the structure once to escape from a '
                        'possible (hard) local minimum.'
                    )
                else:
                    self.report_error_handled(
                        calculation, 'bfgs history failure after the structure was already rattled once, aborting...'
                    )
                    return ProcessHandlerReport(True, self.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE)

                self.set_restart_type(RestartType.FROM_CHARGE_DENSITY, calculation.outputs.remote_folder)
                self.report_error_handled(calculation, action)
                return ProcessHandlerReport(True)

            self.ctx.inputs.parameters['IONS']['ion_dynamics'] = 'damp'
            action = 'bfgs history (ionic only) failure: restarting with `damp` dynamics.'

        elif calculation_type == 'vc-relax' and trust_radius_min > 1.0e-4:
            self.ctx.inputs.parameters['IONS']['trust_radius_ini'] = trust_radius_min  # start close
            new_trust_radius_min = trust_radius_min * self.defaults.delta_factor_trust_radius_min
            self.ctx.inputs.parameters['IONS']['trust_radius_min'] = new_trust_radius_min
            action = f'bfgs history (vc-relax) failure: restarting with `trust_radius_min={new_trust_radius_min:.5f}`.'

        elif calculation_type == 'vc-relax':
            self.ctx.inputs.parameters['IONS']['ion_dynamics'] = 'damp'
            self.ctx.inputs.parameters['CELL']['cell_dynamics'] = 'damp-w'
            action = 'bfgs history (vc-relax) failure: restarting with `damp(-w)` dynamics.'

        else:
            return ProcessHandlerReport(False)

        self.ctx.inputs.structure = calculation.outputs.output_structure

        self.set_restart_type(RestartType.FROM_CHARGE_DENSITY, calculation.outputs.remote_folder)
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True)

    @process_handler(
        priority=560,
        exit_codes=[
            PwCalculation.exit_codes.ERROR_IONIC_CONVERGENCE_NOT_REACHED,
            PwCalculation.exit_codes.ERROR_IONIC_CYCLE_EXCEEDED_NSTEP,
            PwCalculation.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE,
            PwCalculation.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE,
        ],
    )
    def handle_relax_recoverable_ionic_convergence_error(self, calculation):
        """Handle various exit codes for recoverable `relax` calculations with failed ionic convergence.

        These exit codes signify that the ionic convergence thresholds were not met, but the output structure is usable,
        so the solution is to simply restart from scratch but from the output structure.
        """
        self.ctx.inputs.structure = calculation.outputs.output_structure
        action = 'no ionic convergence but clean shutdown: restarting from scratch but using output structure.'

        self.set_restart_type(RestartType.FROM_CHARGE_DENSITY, calculation.outputs.remote_folder)
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True)

    @process_handler(
        priority=555,
        exit_codes=[
            PwCalculation.exit_codes.ERROR_RADIAL_FFT_SIGNIFICANT_VOLUME_CONTRACTION,
        ],
    )
    def handle_vcrelax_recoverable_fft_significant_volume_contraction_error(self, calculation):
        """Handle exit code for recoverable `vc-relax` calculations with significant volume contraction.

        This exit code appears when a cell relaxation produces a significant volume scaling (contraction or expansion).
        This means the pseudopotentials tables must be recalculated. This parameter is controlled by `CELL.cell_factor`.
        The solution, as suggested by the QuantumESPRESSO error itself, is to restart with an increased `cell_factor`.
        We then start from scratch using the last output structure and we double the cell factor.
        """
        self.ctx.inputs.structure = calculation.outputs.output_structure
        self.ctx.inputs.parameters.setdefault('CELL', {})  # as it is not compulsory for ``vc-relax`` calculations
        cell_factor = 2 * self.ctx.inputs.parameters['CELL'].get('cell_factor', 2)
        self.ctx.inputs.parameters['CELL']['cell_factor'] = cell_factor

        self.set_restart_type(RestartType.FROM_SCRATCH)
        action = (
            'significant volume scaling but clean shutdown: '
            f'restarting from scratch using output structure and ``cell_factor = {cell_factor}``.'
        )
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True)

    @process_handler(
        priority=550,
        exit_codes=[
            PwCalculation.exit_codes.ERROR_IONIC_CYCLE_ELECTRONIC_CONVERGENCE_NOT_REACHED,
            PwCalculation.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_FINAL_SCF_FAILED,
        ],
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

    @process_handler(
        priority=410,
        exit_codes=[
            PwCalculation.exit_codes.ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED,
        ],
    )
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

    @process_handler(
        priority=420,
        exit_codes=[
            PwCalculation.exit_codes.WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED,
        ],
    )
    def handle_electronic_convergence_warning(self, calculation):
        """Handle `WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED`: consider finished."""
        self.ctx.is_finished = True
        action = 'electronic convergence not reached but inputs say this is ok: consider finished.'
        self.report_error_handled(calculation, action)
        self.results()  # Call the results method to attach the output nodes
        return ProcessHandlerReport(True, self.exit_codes.WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED)
