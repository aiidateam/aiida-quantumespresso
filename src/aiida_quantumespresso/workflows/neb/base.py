# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO neb.x calculation with automated error handling and restarts."""
from aiida import orm
from aiida.common import AttributeDict, exceptions
from aiida.common.lang import type_check
from aiida.engine import BaseRestartWorkChain, ExitCode, ProcessHandlerReport, process_handler, while_
from aiida.plugins import CalculationFactory, GroupFactory

from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import create_kpoints_from_distance
from aiida_quantumespresso.common.types import ElectronicType, RestartType, SpinType
from aiida_quantumespresso.utils.defaults.calculation import pw as qe_defaults

from ..protocols.utils import ProtocolMixin

NebCalculation = CalculationFactory('quantumespresso.neb')
SsspFamily = GroupFactory('pseudo.family.sssp')
PseudoDojoFamily = GroupFactory('pseudo.family.pseudo_dojo')
CutoffsPseudoPotentialFamily = GroupFactory('pseudo.family.cutoffs')


class NebBaseWorkChain(ProtocolMixin, BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO neb.x calculation with automated error handling and restarts."""

    # pylint: disable=too-many-public-methods, too-many-statements

    _process_class = NebCalculation

    defaults = AttributeDict({
        'qe': qe_defaults,
        'delta_threshold_degauss': 30,
        'delta_factor_degauss': 0.1,
        'delta_factor_mixing_beta': 0.8,
        'delta_factor_max_seconds': 0.95,
        'delta_factor_nbnd': 0.05,
        'delta_minimum_nbnd': 4,
        'delta_factor_trust_radius_min': 0.1,
    })

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.expose_inputs(NebCalculation, namespace='neb', exclude=('pw.kpoints',))
        spec.input('kpoints', valid_type=orm.KpointsData, required=False,
            help='An explicit k-points list or mesh. Either this or `kpoints_distance` has to be provided.')
        spec.input('kpoints_distance', valid_type=orm.Float, required=False,
            help='The minimum desired distance in 1/Å between k-points in reciprocal space. The explicit k-points will '
                 'be generated automatically by a calculation function based on the input structure.')
        spec.input('kpoints_force_parity', valid_type=orm.Bool, required=False,
            help='Optional input when constructing the k-points based on a desired `kpoints_distance`. Setting this to '
                 '`True` will force the k-point mesh to have an even number of points along each lattice vector except '
                 'for any non-periodic directions.')

        spec.outline(
            cls.setup,
            cls.validate_kpoints,
            while_(cls.should_run_process)(
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
        )

        spec.expose_outputs(NebCalculation)

        spec.exit_code(202, 'ERROR_INVALID_INPUT_KPOINTS',
            message='Neither the `kpoints` nor the `kpoints_distance` input was specified.')

        spec.exit_code(300, 'ERROR_UNRECOVERABLE_FAILURE',
            message='The calculation failed with an unidentified unrecoverable error.')
        spec.exit_code(310, 'ERROR_KNOWN_UNRECOVERABLE_FAILURE',
            message='The calculation failed with a known unrecoverable error.')
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
        from aiida_quantumespresso.workflows.protocols.utils import get_magnetization, recursive_merge

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

        if spin_type is SpinType.SPIN_ORBIT and overrides is not None and 'pseudo_family' not in overrides:
            pseudo_family = 'PseudoDojo/0.4/PBEsol/FR/standard/upf'

        natoms = len(structure.sites)

        # Update the parameters based on the protocol inputs
        parameters = inputs['pw']['parameters']

        if overrides and 'pseudos' in overrides.get('pw', {}):

            pseudos = overrides['pw']['pseudos']

            if sorted(pseudos.keys()) != sorted(structure.get_kind_names()):
                raise ValueError(f'`pseudos` override needs one value for each of the {len(structure.kinds)} kinds.')

            system_overrides = overrides['pw'].get('parameters', {}).get('SYSTEM', {})

            if not all(key in system_overrides for key in ('ecutwfc', 'ecutrho')):
                raise ValueError(
                    'When overriding the pseudo potentials, both `ecutwfc` and `ecutrho` cutoffs should be '
                    f'provided in the `overrides`: {overrides}'
                )

        else:
            try:
                pseudo_set = (PseudoDojoFamily, SsspFamily, CutoffsPseudoPotentialFamily)
                pseudo_family = orm.QueryBuilder().append(pseudo_set, filters={'label': pseudo_family}).one()[0]
            except exceptions.NotExistent as exception:
                raise ValueError(
                    f'required pseudo family `{pseudo_family}` is not installed. Please use `aiida-pseudo install` to'
                    'install it.'
                ) from exception

            try:
                parameters['SYSTEM']['ecutwfc'], parameters['SYSTEM'][
                    'ecutrho'] = pseudo_family.get_recommended_cutoffs(structure=structure, unit='Ry')
                pseudos = pseudo_family.get_pseudos(structure=structure)
            except ValueError as exception:
                raise ValueError(
                    f'failed to obtain recommended cutoffs for pseudo family `{pseudo_family}`: {exception}'
                ) from exception

        parameters['CONTROL']['etot_conv_thr'] = natoms * meta_parameters['etot_conv_thr_per_atom']
        parameters['ELECTRONS']['conv_thr'] = natoms * meta_parameters['conv_thr_per_atom']

        # If the structure is 2D periodic in the x-y plane, we set assume_isolate to `2D`
        if structure.pbc == (True, True, False):
            parameters['SYSTEM']['assume_isolated'] = '2D'

        if electronic_type is ElectronicType.INSULATOR:
            parameters['SYSTEM']['occupations'] = 'fixed'
            parameters['SYSTEM'].pop('degauss')
            parameters['SYSTEM'].pop('smearing')

        magnetization = get_magnetization(
            structure=structure,
            z_valences={kind.name: pseudos[kind.name].z_valence for kind in structure.kinds},
            initial_magnetic_moments=initial_magnetic_moments,
            spin_type=spin_type,
        )
        if spin_type is SpinType.COLLINEAR:
            parameters['SYSTEM']['starting_magnetization'] = magnetization['starting_magnetization']
            parameters['SYSTEM']['nspin'] = 2

        if spin_type in [SpinType.SPIN_ORBIT, SpinType.NON_COLLINEAR]:
            parameters['SYSTEM']['starting_magnetization'] = magnetization['starting_magnetization']
            parameters['SYSTEM']['angle1'] = magnetization['angle1']
            parameters['SYSTEM']['angle2'] = magnetization['angle2']
            parameters['SYSTEM']['noncolin'] = True
            parameters['SYSTEM']['nspin'] = 4
            if spin_type == SpinType.SPIN_ORBIT:
                parameters['SYSTEM']['lspinorb'] = True

        # If overrides are provided, they are considered absolute
        if overrides:
            parameter_overrides = overrides.get('pw', {}).get('parameters', {})
            parameters = recursive_merge(parameters, parameter_overrides)

            # if tot_magnetization in overrides , remove starting_magnetization from parameters
            if parameters.get('SYSTEM', {}).get('tot_magnetization') is not None:
                parameters.setdefault('SYSTEM', {}).pop('starting_magnetization', None)

        # pylint: disable=no-member
        builder = cls.get_builder()
        builder.neb['code'] = code
        builder.neb.pw['pseudos'] = pseudos
        # builder.pw['structure'] = structure
        builder.neb.pw['parameters'] = orm.Dict(parameters)

        if 'kpoints' in inputs:
            builder.kpoints = inputs['kpoints']
        else:
            builder.kpoints_distance = orm.Float(inputs['kpoints_distance'])
        builder.kpoints_force_parity = orm.Bool(inputs['kpoints_force_parity'])
        builder.max_iterations = orm.Int(inputs['max_iterations'])
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
        self.ctx.inputs = AttributeDict(self.exposed_inputs(NebCalculation, 'neb'))

        self.ctx.inputs.parameters = self.ctx.inputs.parameters.get_dict()
        self.ctx.inputs.pw.parameters = self.ctx.inputs.pw.parameters.get_dict()
        self.ctx.inputs.pw.parameters.setdefault('CONTROL', {})
        self.ctx.inputs.pw.parameters.setdefault('ELECTRONS', {})
        self.ctx.inputs.pw.parameters.setdefault('SYSTEM', {})

        self.ctx.inputs.pw.parameters['CONTROL'].setdefault('calculation', 'scf')

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
                'structure': self.inputs.neb.last_structure,
                'distance': self.inputs.kpoints_distance,
                'force_parity': self.inputs.get('kpoints_force_parity', orm.Bool(False)),
                'metadata': {
                    'call_link_label': 'create_kpoints_from_distance'
                }
            }
            kpoints = create_kpoints_from_distance(**inputs)  # pylint: disable=unexpected-keyword-arg

        self.ctx.inputs.pw.kpoints = kpoints

    def set_restart_type(self, restart_type, parent_folder=None):
        """Set the restart type for the next iteration."""

        if parent_folder is None and restart_type != RestartType.FROM_SCRATCH:
            raise ValueError('When not restarting from scratch, a `parent_folder` must be provided.')

        if restart_type == RestartType.FROM_SCRATCH:
            self.ctx.inputs.parameters['PATH']['restart_mode'] = 'from_scratch'
            self.ctx.inputs.pop('parent_folder', None)

        elif restart_type == RestartType.FULL:
            self.ctx.inputs.parameters['PATH']['restart_mode'] = 'restart'
            self.ctx.inputs.parent_folder = parent_folder

    def report_error_handled(self, calculation, action):
        """Report an action taken for a calculation that has failed.

        This should be called in a registered error handler if its condition is met and an action was taken.

        :param calculation: the failed calculation node
        :param action: a string message with the action taken
        """
        arguments = [calculation.process_label, calculation.pk, calculation.exit_status, calculation.exit_message]
        self.report('{}<{}> failed with exit status {}: {}'.format(*arguments))
        self.report(f'Action taken: {action}')

    @process_handler(
        priority=585,
        exit_codes=[
            NebCalculation.exit_codes.ERROR_COMPUTING_CHOLESKY,
            NebCalculation.exit_codes.ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED,
            NebCalculation.exit_codes.ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE,
            NebCalculation.exit_codes.ERROR_ZHEGVD_FAILED,
            NebCalculation.exit_codes.ERROR_QR_FAILED,
            NebCalculation.exit_codes.ERROR_EIGENVECTOR_CONVERGENCE,
            NebCalculation.exit_codes.ERROR_BROYDEN_FACTORIZATION,
        ]
    )
    def handle_diagonalization_errors(self, calculation):
        """Handle known issues related to the diagonalization.

        We use the following strategy. When a diagonalization algorithm fails, we try using an other one
        still not used. Conjugate gradient (CG) is kept as last option, as it is the slowest among the
        available ones, but on the contrary it is the most stable as well, thus kept as `last resort`.

        Once the error handler has tried all ``diagonalization`` options, abort.
        """
        current = self.ctx.inputs.pw.parameters['ELECTRONS'].get('diagonalization', 'david')

        if 'diagonalizations' not in self.ctx:
            # Initialize a list to track diagonalisations that haven't been tried in reverse order or preference
            self.ctx.diagonalizations = [value for value in ['cg', 'paro', 'ppcg', 'david'] if value != current.lower()]

        try:
            new = self.ctx.diagonalizations.pop()
            self.ctx.inputs.pw.parameters['ELECTRONS']['diagonalization'] = new
            action = f'found diagonalization issues for ``{current}``, switching to ``{new}`` diagonalization.'
            self.report_error_handled(calculation, action)
            return ProcessHandlerReport(True)
        except IndexError:
            action = 'found diagonalization issues but already exploited all supported algorithms, aborting...'
            self.report_error_handled(calculation, action)
            return ProcessHandlerReport(True, self.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE)

    @process_handler(priority=575, exit_codes=[
        NebCalculation.exit_codes.ERROR_NEB_INTERRUPTED_PARTIAL_TRAJECTORY,
    ])
    def handle_neb_interrupted_partial_trajectory(self, calculation):
        """Handle `ERROR_NEB_INTERRUPTED_PARTIAL_TRAJECTORY` and exit code.

        In this case the calculation shut down cleanly and we can do a full restart.
        """

        self.set_restart_type(RestartType.FULL, calculation.outputs.remote_folder)
        self.report_error_handled(calculation, "restarting in full with `CONTROL.restart_mode` = 'restart'")

        return ProcessHandlerReport(True)

    @process_handler(priority=575, exit_codes=[
        NebCalculation.exit_codes.ERROR_NEB_CYCLE_EXCEEDED_NSTEP,
    ])
    def handle_neb_cycle_exceeded_nstep(self, calculation):
        """Handle `ERROR_NEB_CYCLE_EXCEEDED_NSTEP` and exit code.

        In this case the calculation shut down cleanly and we can do a full restart.
        """
        self.ctx.inputs.parameters['PATH'].setdefault('nstep_path', 1)
        input_nsteps = self.inputs.neb.parameters['PATH']['nstep_path'] if 'nstep_path' in self.inputs.neb.parameters[
            'PATH'] else 1
        self.ctx.inputs.parameters['PATH']['nstep_path'] += input_nsteps

        self.set_restart_type(RestartType.FULL, calculation.outputs.remote_folder)
        self.report_error_handled(calculation, "restarting in full with `CONTROL.restart_mode` = 'restart'")

        return ProcessHandlerReport(True)

    @process_handler(priority=410, exit_codes=[
        NebCalculation.exit_codes.ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED,
    ])
    def handle_electronic_convergence_not_reached(self, calculation):
        """Handle `ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED` error.

        Decrease the mixing beta and fully restart from the previous calculation.
        """
        factor = self.defaults.delta_factor_mixing_beta
        mixing_beta = self.ctx.inputs.pw.parameters.get('ELECTRONS',
                                                        {}).get('mixing_beta', self.defaults.qe.mixing_beta)
        mixing_beta_new = mixing_beta * factor

        self.ctx.inputs.pw.parameters['ELECTRONS']['mixing_beta'] = mixing_beta_new
        action = f'reduced beta mixing from {mixing_beta} to {mixing_beta_new} and restarting from the last calculation'

        self.set_restart_type(RestartType.FULL, calculation.outputs.remote_folder)
        self.report_error_handled(calculation, action)
        return ProcessHandlerReport(True)

    @process_handler(priority=600)
    def handle_unrecoverable_failure(self, calculation):
        """Handle calculations with an exit status below 400 which are unrecoverable, so abort the work chain."""
        if calculation.is_failed and calculation.exit_status < 400:
            self.report_error_handled(calculation, 'unrecoverable error, aborting...')
            return ProcessHandlerReport(True, self.exit_codes.ERROR_UNRECOVERABLE_FAILURE)
