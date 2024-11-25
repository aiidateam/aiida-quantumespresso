"""Workchain to relax a structure using Quantum ESPRESSO pw.x."""

from aiida import orm
from aiida.common import AttributeDict, exceptions
from aiida.common.lang import type_check
from aiida.engine import ToContext, WorkChain, append_, if_, while_
from aiida.orm import StructureData as LegacyStructureData
from aiida.plugins import DataFactory

from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import create_kpoints_from_distance
from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.common.types import RelaxType
from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

from ..protocols.utils import ProtocolMixin

try:
    StructureData = DataFactory('atomistic.structure')
except exceptions.MissingEntryPointError:
    structures_classes = (LegacyStructureData,)
else:
    structures_classes = (LegacyStructureData, StructureData)


def validate_inputs(inputs, _):
    """Validate the top level namespace."""
    parameters = inputs['base']['pw']['parameters'].get_dict()

    if 'calculation' not in parameters.get('CONTROL', {}):
        return 'The parameters in `base.pw.parameters` do not specify the required key `CONTROL.calculation`.'


class PwRelaxWorkChain(ProtocolMixin, WorkChain):
    """Workchain to relax a structure using Quantum ESPRESSO pw.x."""

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        super().define(spec)
        spec.expose_inputs(
            PwBaseWorkChain,
            namespace='base_init_relax',
            exclude=('clean_workdir', 'pw.structure', 'pw.parent_folder'),
            namespace_options={
                'required': False,
                'populate_defaults': False,
                'help': (
                    'Inputs for the `PwBaseWorkChain` that runs an initial geometry optimization, typically with looser'
                    'precision settings to find a geometry that is already closer to the final one quickly.'
                ),
            },
        )
        spec.expose_inputs(
            PwBaseWorkChain,
            namespace='base_relax',
            exclude=('clean_workdir', 'pw.structure', 'pw.parent_folder'),
            namespace_options={'help': 'Inputs for the `PwBaseWorkChain` for the main relax loop.'})
        spec.input('structure', valid_type=structures_classes, help='The inputs structure.')
        spec.input('meta_convergence', valid_type=orm.Bool, default=lambda: orm.Bool(True),
            help='If `True` the workchain will perform a meta-convergence on the cell volume.')
        spec.input('max_meta_convergence_iterations', valid_type=orm.Int, default=lambda: orm.Int(5),
            help='The maximum number of variable cell relax iterations in the meta convergence cycle.')
        spec.input('volume_convergence', valid_type=orm.Float, default=lambda: orm.Float(0.01),
            help='The volume difference threshold between two consecutive meta convergence iterations.')
        spec.input('clean_workdir', valid_type=orm.Bool, default=lambda: orm.Bool(False),
            help='If `True`, work directories of all called calculation will be cleaned at the end of execution.')
        spec.inputs.validator = cls.validate_inputs
        spec.outline(
            cls.setup,
            if_(cls.should_run_init_relax)(
                cls.run_init_relax,
                cls.inspect_init_relax,
            ),
            while_(cls.should_run_relax)(
                cls.run_relax,
                cls.inspect_relax,
            ),
            cls.results,
        )
        spec.exit_code(
            400,
            'ERROR_MAX_ITERATIONS_EXCEEDED',
            message='The maximum number of meta convergence iterations was exceeded.',
        )
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_RELAX', message='the relax `PwBaseWorkChain` sub process failed')
        spec.exit_code(
            402,
            'ERROR_SUB_PROCESS_FAILED_FINAL_SCF',
            message='the final scf PwBaseWorkChain sub process failed '
            '(deprecated since the final SCF has been removed)',
        )
        spec.exit_code(
            403,
            'ERROR_SUB_PROCESS_FAILED_INIT_RELAX',
            message='the initial relaxation `PwBaseWorkChain` sub process failed',
        )
        spec.expose_outputs(PwBaseWorkChain, exclude=('output_structure',))
<<<<<<< HEAD
        spec.output(
            'output_structure', valid_type=orm.StructureData, required=False, help='The successfully relaxed structure.'
        )
=======
        spec.output('output_structure', valid_type=structures_classes, required=False,
            help='The successfully relaxed structure.')
        # yapf: enable
>>>>>>> 76cb491 (Optional support for aiida-atomistic)

    @staticmethod
    def validate_inputs(inputs, _):
        """Validate the top level namespace."""
        parameters = inputs['base_relax']['pw']['parameters'].get_dict()

        if 'calculation' not in parameters.get('CONTROL', {}):
            return 'The parameters in `base.pw.parameters` do not specify the required key `CONTROL.calculation`.'

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        from ..protocols import pw as pw_protocols

        return files(pw_protocols) / 'relax.yaml'

    @classmethod
    def get_builder_from_protocol(
        cls, code, structure, protocol=None, overrides=None, relax_type=RelaxType.POSITIONS_CELL, options=None, **kwargs
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param structure: the ``StructureData`` instance to use.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param relax_type: the relax type to use: should be a value of the enum ``common.types.RelaxType``.
        :param options: A dictionary of options that will be recursively set for the ``metadata.options`` input of all
            the ``CalcJobs`` that are nested in this work chain.
        :param kwargs: additional keyword arguments that will be passed to the ``get_builder_from_protocol`` of all the
            sub processes that are called by this workchain.
        :return: a process builder instance with all inputs defined ready for launch.
        """
        type_check(relax_type, RelaxType)

        inputs = cls.get_protocol_inputs(protocol, overrides)

        args = (code, structure, protocol)
        base_relax = PwBaseWorkChain.get_builder_from_protocol(
            *args, overrides=inputs.get('base_relax', None), options=options, **kwargs
        )
        base_init_relax = PwBaseWorkChain.get_builder_from_protocol(
            *args, overrides=inputs.get('base_init_relax', None), options=options, **kwargs
        )
        base_relax['pw'].pop('structure', None)
        base_relax.pop('clean_workdir', None)
        base_init_relax['pw'].pop('structure', None)
        base_init_relax.pop('clean_workdir', None)

        for namespace in (base_relax, base_init_relax):
            # Quantum ESPRESSO currently only supports optimization of the volume for simple cubic systems. It requires
            # to set `ibrav=1` or the code will except.
            if relax_type in (RelaxType.VOLUME, RelaxType.POSITIONS_VOLUME):
                raise ValueError(f'relax type `{relax_type} is not yet supported.')

            if relax_type in (RelaxType.VOLUME, RelaxType.SHAPE, RelaxType.CELL):
                namespace.pw.settings = orm.Dict(
                    PwRelaxWorkChain._fix_atomic_positions(structure, base_relax.pw.settings)
                )

            if relax_type is RelaxType.NONE:
                namespace.pw.parameters['CONTROL']['calculation'] = 'scf'
                namespace.pw.parameters.base.attributes.delete('CELL')

            elif relax_type is RelaxType.POSITIONS:
                namespace.pw.parameters['CONTROL']['calculation'] = 'relax'
                namespace.pw.parameters.base.attributes.delete('CELL')
            else:
                namespace.pw.parameters['CONTROL']['calculation'] = 'vc-relax'

            if relax_type in (RelaxType.VOLUME, RelaxType.POSITIONS_VOLUME):
                namespace.pw.parameters['CELL']['cell_dofree'] = 'volume'

            if relax_type in (RelaxType.SHAPE, RelaxType.POSITIONS_SHAPE):
                namespace.pw.parameters['CELL']['cell_dofree'] = 'shape'

            if relax_type in (RelaxType.CELL, RelaxType.POSITIONS_CELL):
                pbc_cell_dofree_map = {
                    (True, True, True): 'all',
                    (True, False, False): 'x',
                    (False, True, False): 'y',
                    (False, False, True): 'z',
                    (True, True, False): '2Dxy',
                }
                if structure.pbc in pbc_cell_dofree_map:
                    namespace.pw.parameters['CELL']['cell_dofree'] = pbc_cell_dofree_map[structure.pbc]
                else:
                    raise ValueError(
                        f'Structures with periodic boundary conditions `{structure.pbc}` are not supported.'
                    )

        builder = cls.get_builder()
        builder.base_relax = base_relax
        builder.base_init_relax = base_init_relax
        builder.structure = structure
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])
        builder.max_meta_convergence_iterations = orm.Int(inputs['max_meta_convergence_iterations'])
        builder.meta_convergence = orm.Bool(inputs['meta_convergence'])

        return builder

    def setup(self):
        """Input validation and context setup."""
        self.ctx.current_structure = self.inputs.structure
        self.ctx.current_number_of_bands = None
        self.ctx.iteration = 0

        self.ctx.relax_inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='base_relax'))
        self.ctx.relax_inputs.pw.parameters = self.ctx.relax_inputs.pw.parameters.get_dict()

        self.ctx.relax_inputs.pw.parameters.setdefault('CONTROL', {})

        # Set the meta_convergence and add it to the context
        self.ctx.meta_convergence = self.inputs.meta_convergence.value
        volume_cannot_change = (
            self.ctx.relax_inputs.pw.parameters['CONTROL'].get('calculation', 'scf') in ('scf', 'relax')
            or self.ctx.relax_inputs.pw.parameters.get('CELL', {}).get('cell_dofree', None) == 'shape'
        )
        if self.ctx.meta_convergence and volume_cannot_change:
            self.report(
                'No change in volume possible for the provided base input parameters. Meta convergence is turned off.'
            )
            self.ctx.meta_convergence = False
        # Already converged if volume can not change
        self.ctx.converged = volume_cannot_change

    def should_run_init_relax(self):
        """Return whether an initial relaxation should be run."""
        return 'base_init_relax' in self.inputs

    def run_init_relax(self):
        """Run the `PwBaseWorkChain` to run an initial relaxation."""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='base_init_relax'))
        inputs.pw.structure = self.ctx.current_structure

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        inputs.metadata.call_link_label = 'init_relax'

        base_wc = self.submit(PwBaseWorkChain, **inputs)
        self.report(f'launching PwBaseWorkChain<{base_wc.pk}> for initial relaxation.')

        return ToContext(base_init_relax_workchain=base_wc)

    def inspect_init_relax(self):
        """Inspect the result of the initial relax `PwBaseWorkChain`."""
        workchain = self.ctx.base_init_relax_workchain

        if not workchain.is_finished_ok:
            self.report(f'final scf PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_INIT_RELAX

        self.ctx.current_structure = workchain.outputs.output_structure
        self.ctx.current_number_of_bands = workchain.outputs.output_parameters.get_dict()['number_of_bands']

    def should_run_relax(self):
        """Return whether a relaxation workchain should be run.

        This is the case as long as either:

        1. There are still Pulay stresses in the system. When running a `vc-relax`, Quantum ESPRESSO will automatically
           run a final SCF after the geometry optimization, where it regenerates the plane wave basis sets using the
           final lattice of the final geometry. If the stresses are too large in this final SCF, this will means there
           are still Pulay stresses and the `PwBaseWorkChain` will check for this and return exit code 501.
        2. The k-points mesh was defined as a density using `kpoints_distance` for the `base_relax` input name space,
           and the unit cell has changed in such a way that regenerating the k-points mesh using this density increases
           the mesh compared to the previous run.

        If one of these two conditions is met, another relaxation should be run as long as the maximum number of meta
        convergence iterations is not exceeded.
        """
        # If no work chain has been run, we should at least run one
        if 'base_relax_workchains' not in self.ctx:
            return True

        # If meta convergence is switched off we are done
        if not self.ctx.meta_convergence:
            self.report('Meta-convergence disabled: workchain completed after a single iteration.')
            return False

        # Stop if the maximum number of meta iterations has been reached
        if self.ctx.iteration == self.inputs.max_meta_convergence_iterations.value:
            return False

        base_relax_workchain = self.ctx.base_relax_workchains[-1]

        # If the last work chain still found Pulay stresses in the final SCF, continue
        pulay_exit_status = PwCalculation.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF.status
        if base_relax_workchain.exit_status == pulay_exit_status:
            self.report('Pulay stresses still present, running another geometry optimization.')
            return True

        # If the kpoints are defined as a density, make sure the kpoints mesh is the same for the new structur
        if 'kpoints_distance' in self.ctx.relax_inputs:
            input_kpts_mesh, _ = (
                base_relax_workchain.base.links.get_outgoing(link_label_filter='iteration_01')
                .first()
                .node.inputs.kpoints.get_kpoints_mesh()
            )
            inputs_create_kpoints = {
                'structure': base_relax_workchain.outputs.output_structure,
                'distance': self.ctx.relax_inputs.kpoints_distance,
                'force_parity': self.ctx.relax_inputs.get('kpoints_force_parity', orm.Bool(False)),
                'metadata': {'store_provenance': False},
            }
            new_kpts_mesh, _ = create_kpoints_from_distance(**inputs_create_kpoints).get_kpoints_mesh()

            if not all(k1 <= k2 for k1, k2 in zip(new_kpts_mesh, input_kpts_mesh)):
                self.report(
                    'Density of k-points mesh has increased for the specified `kpoints_distance` due to a change in the'
                    ' unit cell. Running another geometry optimization with new mesh.'
                )
                return True

        self.report(f'Work chain completed after {self.ctx.iteration} iterations.')
        self.ctx.converged = True
        return False

    def run_relax(self):
        """Run the `PwBaseWorkChain` to run a relax `PwCalculation`."""
        self.ctx.iteration += 1
        inputs = self.ctx.relax_inputs
        inputs.pw.structure = self.ctx.current_structure

        # If one of the nested `PwBaseWorkChain`s changed the number of bands, apply it here
        if self.ctx.current_number_of_bands is not None:
            inputs.pw.parameters.setdefault('SYSTEM', {})['nbnd'] = self.ctx.current_number_of_bands

        # Set the `CALL` link label
        inputs.metadata.call_link_label = f'iteration_{self.ctx.iteration:02d}'
        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)

        base_wc = self.submit(PwBaseWorkChain, **inputs)
        self.report(f'launching PwBaseWorkChain<{base_wc.pk}>')

        return ToContext(base_relax_workchains=append_(base_wc))

    def inspect_relax(self):
        """Inspect the results of the last `PwBaseWorkChain`."""
        workchain = self.ctx.base_relax_workchains[-1]

        if workchain.is_excepted or workchain.is_killed:
            self.report('relax PwBaseWorkChain was excepted or killed')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        # The following list of `PwBaseWorkChain` exit status should not interrupt the work chain
        acceptable_statuses = ['ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF']

        if workchain.is_failed and workchain.exit_status not in PwBaseWorkChain.get_exit_statuses(acceptable_statuses):
            self.report(f'relax PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        try:
            structure = workchain.outputs.output_structure
        except exceptions.NotExistent:
            # If the calculation is set to 'scf', this is expected, so we are done
            if self.ctx.relax_inputs.pw.parameters['CONTROL']['calculation'] == 'scf':
                return

            self.report('`vc-relax` or `relax` PwBaseWorkChain finished successfully but without output structure')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        # Set relaxed structure as input structure for next iteration
        self.ctx.current_structure = structure
        self.ctx.current_number_of_bands = workchain.outputs.output_parameters.get_dict()['number_of_bands']

    def results(self):
        """Attach the output parameters and structure of the last workchain to the outputs."""
        if self.ctx.iteration == self.inputs.max_meta_convergence_iterations.value and not self.ctx.converged:
            self.report('Maximum number of meta convergence iterations reached.')
            return self.exit_codes.ERROR_MAX_ITERATIONS_EXCEEDED

        # Get the latest relax workchain and pass the outputs
        final_relax_workchain = self.ctx.base_relax_workchains[-1]

        if self.ctx.relax_inputs.pw.parameters['CONTROL']['calculation'] != 'scf':
            self.out('output_structure', final_relax_workchain.outputs.output_structure)

        self.out_many(self.exposed_outputs(final_relax_workchain, PwBaseWorkChain))

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super().on_terminated()

        if self.inputs.clean_workdir.value is False:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = []

        for called_descendant in self.node.called_descendants:
            if isinstance(called_descendant, orm.CalcJobNode):
                try:
                    called_descendant.outputs.remote_folder._clean()  # noqa: SLF001
                    cleaned_calcs.append(called_descendant.pk)
                except (OSError, KeyError):
                    pass

        if cleaned_calcs:
            self.report(f'cleaned remote folders of calculations: {" ".join(map(str, cleaned_calcs))}')

    @staticmethod
    def _fix_atomic_positions(structure, settings):
        """Fix the atomic positions, by setting the `FIXED_COORDS` key in the `settings` input node."""
        settings = settings.get_dict() if settings is not None else {}

        settings['FIXED_COORDS'] = [[True, True, True]] * len(structure.sites)

        return settings
