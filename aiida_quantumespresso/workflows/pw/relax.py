# -*- coding: utf-8 -*-
"""Workchain to relax a structure using Quantum ESPRESSO pw.x."""
from aiida import orm
from aiida.common import AttributeDict, exceptions
from aiida.engine import WorkChain, ToContext, if_, while_, append_
from aiida.plugins import CalculationFactory, WorkflowFactory

from aiida_quantumespresso.common.types import RelaxType
from aiida_quantumespresso.utils.mapping import prepare_process_inputs

from ..protocols.utils import ProtocolMixin

PwCalculation = CalculationFactory('quantumespresso.pw')
PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')


def validate_final_scf(value, _):
    """Validate the final scf input."""
    if isinstance(value, orm.Bool) and value:
        import warnings
        from aiida.common.warnings import AiidaDeprecationWarning
        warnings.warn(
            'this input is deprecated and will be removed. If you want to run a final scf, specify the inputs that '
            'should be used in the `base_final_scf` namespace.', AiidaDeprecationWarning
        )


def validate_relaxation_scheme(value, _):
    """Validate the relaxation scheme input."""
    if value:
        import warnings
        from aiida.common.warnings import AiidaDeprecationWarning
        warnings.warn(
            'the `relaxation_scheme` input is deprecated and will be removed. Use the ``relax_type`` input instead. '
            'Accepted values are values of the ``aiida_quantumespresso.common.types.RelaxType`` enum',
            AiidaDeprecationWarning
        )


def validate_relax_type(value, _):
    """Validate the relax type input."""
    if value:
        try:
            RelaxType(value.value)
        except ValueError:
            return f'`{value.value}` is not a valid value of `RelaxType`.'


class PwRelaxWorkChain(ProtocolMixin, WorkChain):
    """Workchain to relax a structure using Quantum ESPRESSO pw.x."""

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.expose_inputs(PwBaseWorkChain, namespace='base',
            exclude=('clean_workdir', 'pw.structure', 'pw.parent_folder'),
            namespace_options={'help': 'Inputs for the `PwBaseWorkChain` for the main relax loop.'})
        spec.expose_inputs(PwBaseWorkChain, namespace='base_final_scf',
            exclude=('clean_workdir', 'pw.structure', 'pw.parent_folder'),
            namespace_options={'required': False, 'populate_defaults': False,
                'help': 'Inputs for the `PwBaseWorkChain` for the final scf.'})
        spec.input('structure', valid_type=orm.StructureData, help='The inputs structure.')
        spec.input('final_scf', valid_type=orm.Bool, default=lambda: orm.Bool(False), validator=validate_final_scf,
            help='If `True`, a final SCF calculation will be performed on the successfully relaxed structure.')
        spec.input('relaxation_scheme', valid_type=orm.Str, required=False, validator=validate_relaxation_scheme,
            help='The relaxation scheme to use: choose either `relax` or `vc-relax` for variable cell relax.')
        spec.input('relax_type', valid_type=orm.Str, default=lambda: orm.Str(RelaxType.ATOMS_CELL.value),
            validator=validate_relax_type,
            help='The relax type to use: should be a value of the enum ``common.types.RelaxType``.')
        spec.input('meta_convergence', valid_type=orm.Bool, default=lambda: orm.Bool(True),
            help='If `True` the workchain will perform a meta-convergence on the cell volume.')
        spec.input('max_meta_convergence_iterations', valid_type=orm.Int, default=lambda: orm.Int(5),
            help='The maximum number of variable cell relax iterations in the meta convergence cycle.')
        spec.input('volume_convergence', valid_type=orm.Float, default=lambda: orm.Float(0.01),
            help='The volume difference threshold between two consecutive meta convergence iterations.')
        spec.input('clean_workdir', valid_type=orm.Bool, default=lambda: orm.Bool(False),
            help='If `True`, work directories of all called calculation will be cleaned at the end of execution.')
        spec.outline(
            cls.setup,
            while_(cls.should_run_relax)(
                cls.run_relax,
                cls.inspect_relax,
            ),
            if_(cls.should_run_final_scf)(
                cls.run_final_scf,
                cls.inspect_final_scf,
            ),
            cls.results,
        )
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_RELAX',
            message='the relax PwBaseWorkChain sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_FINAL_SCF',
            message='the final scf PwBaseWorkChain sub process failed')
        spec.expose_outputs(PwBaseWorkChain, exclude=('output_structure',))
        spec.output('output_structure', valid_type=orm.StructureData, required=False,
            help='The successfully relaxed structure, unless `relax_type is RelaxType.NONE`.')
        # yapf: enable

    @classmethod
    def get_builder_from_protocol(cls, code, structure, protocol=None, overrides=None, **kwargs):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param structure: the ``StructureData`` instance to use.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param kwargs: additional keyword arguments that will be passed to the ``get_builder_from_protocol`` of all the
            sub processes that are called by this workchain.
        :return: a process builder instance with all inputs defined ready for launch.
        """
        args = (code, structure, protocol)
        inputs = cls.get_protocol_inputs(protocol, overrides)
        builder = cls.get_builder()

        base = PwBaseWorkChain.get_builder_from_protocol(*args, overrides=inputs.get('base', None), **kwargs)
        base_final_scf = PwBaseWorkChain.get_builder_from_protocol(
            *args, overrides=inputs.get('base_final_scf', None), **kwargs
        )

        base['pw'].pop('structure', None)
        base.pop('clean_workdir', None)
        base_final_scf['pw'].pop('structure', None)
        base_final_scf.pop('clean_workdir', None)

        builder.base = base
        builder.base_final_scf = base_final_scf
        builder.structure = structure
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])
        builder.max_meta_convergence_iterations = orm.Int(inputs['max_meta_convergence_iterations'])
        builder.meta_convergence = orm.Bool(inputs['meta_convergence'])
        builder.relax_type = orm.Str(inputs['relax_type'])
        builder.volume_convergence = orm.Float(inputs['volume_convergence'])

        return builder

    def setup(self):
        """Input validation and context setup."""
        self.ctx.current_number_of_bands = None
        self.ctx.current_structure = self.inputs.structure
        self.ctx.current_cell_volume = None
        self.ctx.is_converged = False
        self.ctx.iteration = 0

        # Determine the correct `RelaxType` and add it to the context
        if 'relaxation_scheme' in self.inputs:
            if self.inputs.relaxation_scheme.value == 'relax':
                self.ctx.relax_type = RelaxType.ATOMS
            elif self.inputs.relaxation_scheme.value == 'vc-relax':
                self.ctx.relax_type = RelaxType.ATOMS_CELL
            else:
                raise ValueError('unsupported value for the `relaxation_scheme` input.')
        else:
            self.ctx.relax_type = RelaxType(self.inputs.relax_type)

        # Set the meta_convergence and add it to the context
        self.ctx.meta_convergence = self.inputs.meta_convergence.value
        if self.inputs.meta_convergence and self.ctx.relax_type in [RelaxType.NONE, RelaxType.ATOMS, RelaxType.SHAPE]:
            self.report(
                f'Since there is no change in volume for `{self.ctx.relax_type}`, meta convergence is turned off.'
            )
            self.ctx.meta_convergence = False

        # Add the final scf inputs to the context if a final scf should be run
        if self.inputs.final_scf and 'base_final_scf' in self.inputs:
            raise ValueError('cannot specify `final_scf=True` and `base_final_scf` at the same time.')
        elif self.inputs.final_scf:
            self.ctx.final_scf_inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='base'))
        elif 'base_final_scf' in self.inputs:
            self.ctx.final_scf_inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='base_final_scf'))

        if self.ctx.relax_type is RelaxType.NONE:
            self.report('Work chain will not run final SCF for `RelaxType.NONE`.')
            self.ctx.pop('final_scf_inputs')

    def should_run_relax(self):
        """Return whether a relaxation workchain should be run.

        This is the case as long as the volume change between two consecutive relaxation runs is larger than the volume
        convergence threshold value and the maximum number of meta convergence iterations is not exceeded.
        """
        return not self.ctx.is_converged and self.ctx.iteration < self.inputs.max_meta_convergence_iterations.value

    def should_run_final_scf(self):
        """Return whether after successful relaxation a final scf calculation should be run.

        If the maximum number of meta convergence iterations has been exceeded and convergence has not been reached, the
        structure cannot be considered to be relaxed and the final scf should not be run.
        """
        return self.ctx.is_converged and 'final_scf_inputs' in self.ctx

    def run_relax(self):
        """Run the `PwBaseWorkChain` to run a relax `PwCalculation`."""
        self.ctx.iteration += 1

        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='base'))
        inputs.pw.structure = self.ctx.current_structure
        inputs.pw.parameters = inputs.pw.parameters.get_dict()

        inputs.pw.parameters.setdefault('CELL', {})
        inputs.pw.parameters.setdefault('CONTROL', {})
        inputs.pw.parameters['CONTROL']['restart_mode'] = 'from_scratch'

        if self.ctx.relax_type in [RelaxType.VOLUME, RelaxType.SHAPE, RelaxType.CELL]:
            inputs.pw.settings = self._fix_atomic_positions(inputs.pw.structure, inputs.pw.get('settings', None))

        if self.ctx.relax_type is RelaxType.NONE:
            inputs.pw.parameters['CONTROL']['calculation'] = 'scf'
            inputs.pw.parameters.pop('CELL', None)
        elif self.ctx.relax_type is RelaxType.ATOMS:
            inputs.pw.parameters['CONTROL']['calculation'] = 'relax'
            inputs.pw.parameters.pop('CELL', None)
        else:
            inputs.pw.parameters['CONTROL']['calculation'] = 'vc-relax'

        if self.ctx.relax_type in [RelaxType.VOLUME, RelaxType.ATOMS_VOLUME]:
            inputs.pw.parameters['CELL']['cell_dofree'] = 'volume'

        if self.ctx.relax_type in [RelaxType.SHAPE, RelaxType.ATOMS_SHAPE]:
            inputs.pw.parameters['CELL']['cell_dofree'] = 'shape'

        if self.ctx.relax_type in [RelaxType.CELL, RelaxType.ATOMS_CELL]:
            inputs.pw.parameters['CELL']['cell_dofree'] = 'all'

        # If one of the nested `PwBaseWorkChains` changed the number of bands, apply it here
        if self.ctx.current_number_of_bands is not None:
            inputs.pw.parameters.setdefault('SYSTEM', {})['nbnd'] = self.ctx.current_number_of_bands

        # Set the `CALL` link label
        inputs.metadata.call_link_label = f'iteration_{self.ctx.iteration:02d}'

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launching PwBaseWorkChain<{running.pk}>')

        return ToContext(workchains=append_(running))

    def inspect_relax(self):
        """Inspect the results of the last `PwBaseWorkChain`.

        Compare the cell volume of the relaxed structure of the last completed workchain with the previous. If the
        difference ratio is less than the volume convergence threshold we consider the cell relaxation converged.
        """
        workchain = self.ctx.workchains[-1]

        acceptable_statuses = ['ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF']

        if workchain.is_excepted or workchain.is_killed:
            self.report('relax PwBaseWorkChain was excepted or killed')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        if workchain.is_failed and workchain.exit_status not in PwBaseWorkChain.get_exit_statuses(acceptable_statuses):
            self.report(f'relax PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        try:
            structure = workchain.outputs.output_structure
        except exceptions.NotExistent:
            # If relax_type is RelaxType.NONE, this is expected, so we are done
            if self.ctx.relax_type is RelaxType.NONE:
                self.ctx.is_converged = True
                return

            self.report('`vc-relax` or `relax` PwBaseWorkChain finished successfully but without output structure')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        prev_cell_volume = self.ctx.current_cell_volume
        curr_cell_volume = structure.get_cell_volume()

        # Set relaxed structure as input structure for next iteration
        self.ctx.current_structure = structure
        self.ctx.current_number_of_bands = workchain.outputs.output_parameters.get_dict()['number_of_bands']
        self.report(f'after iteration {self.ctx.iteration} cell volume of relaxed structure is {curr_cell_volume}')

        # After first iteration, simply set the cell volume and restart the next base workchain
        if not prev_cell_volume:
            self.ctx.current_cell_volume = curr_cell_volume

            # If meta convergence is switched off we are done
            if not self.ctx.meta_convergence:
                self.ctx.is_converged = True
            return

        # Check whether the cell volume is converged
        volume_threshold = self.inputs.volume_convergence.value
        volume_difference = abs(prev_cell_volume - curr_cell_volume) / prev_cell_volume

        if volume_difference < volume_threshold:
            self.ctx.is_converged = True
            self.report(
                'relative cell volume difference {} smaller than convergence threshold {}'.format(
                    volume_difference, volume_threshold
                )
            )
        else:
            self.report(
                'current relative cell volume difference {} larger than convergence threshold {}'.format(
                    volume_difference, volume_threshold
                )
            )

        self.ctx.current_cell_volume = curr_cell_volume

        return

    def run_final_scf(self):
        """Run the `PwBaseWorkChain` to run a final scf `PwCalculation` for the relaxed structure."""
        inputs = self.ctx.final_scf_inputs
        inputs.pw.structure = self.ctx.current_structure
        inputs.pw.parameters = inputs.pw.parameters.get_dict()

        inputs.pw.parameters.setdefault('CONTROL', {})
        inputs.pw.parameters['CONTROL']['calculation'] = 'scf'
        inputs.pw.parameters['CONTROL']['restart_mode'] = 'from_scratch'
        inputs.pw.parameters.pop('CELL', None)
        inputs.metadata.call_link_label = 'final_scf'

        if self.ctx.current_number_of_bands is not None:
            inputs.pw.parameters.setdefault('SYSTEM', {})['nbnd'] = self.ctx.current_number_of_bands

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launching PwBaseWorkChain<{running.pk}> for final scf')

        return ToContext(workchain_scf=running)

    def inspect_final_scf(self):
        """Inspect the result of the final scf `PwBaseWorkChain`."""
        workchain = self.ctx.workchain_scf

        if not workchain.is_finished_ok:
            self.report(f'final scf PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_FINAL_SCF

    def results(self):
        """Attach the output parameters and structure of the last workchain to the outputs."""
        if self.ctx.is_converged and self.ctx.iteration <= self.inputs.max_meta_convergence_iterations.value:
            self.report(f'workchain completed after {self.ctx.iteration} iterations')
        else:
            self.report('maximum number of meta convergence iterations exceeded')

        # Get the latest relax workchain and pass the outputs
        final_relax_workchain = self.ctx.workchains[-1]

        if self.ctx.relax_type is not RelaxType.NONE:
            self.out('output_structure', final_relax_workchain.outputs.output_structure)

        if 'final_scf_inputs' in self.ctx:
            self.out_many(self.exposed_outputs(self.ctx.workchain_scf, PwBaseWorkChain))
        else:
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
                    called_descendant.outputs.remote_folder._clean()  # pylint: disable=protected-access
                    cleaned_calcs.append(called_descendant.pk)
                except (IOError, OSError, KeyError):
                    pass

        if cleaned_calcs:
            self.report(f"cleaned remote folders of calculations: {' '.join(map(str, cleaned_calcs))}")

    @staticmethod
    def _fix_atomic_positions(structure, settings):
        """Fix the atomic positions, by setting the `FIXED_COORDS` key in the `settings` input node."""
        if settings is not None:
            settings = settings.get_dict()
        else:
            settings = {}

        settings['FIXED_COORDS'] = [[True, True, True]] * len(structure.sites)

        return settings
