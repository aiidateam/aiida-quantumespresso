# -*- coding: utf-8 -*-
"""Workchain to relax a structure using Quantum ESPRESSO pw.x."""
from __future__ import absolute_import

from six.moves import map

from aiida import orm
from aiida.common import AttributeDict, exceptions
from aiida.engine import WorkChain, ToContext, if_, while_, append_
from aiida.plugins import CalculationFactory, WorkflowFactory
from aiida_quantumespresso.utils.mapping import prepare_process_inputs

PwCalculation = CalculationFactory('quantumespresso.pw')
PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')


class PwRelaxWorkChain(WorkChain):
    """Workchain to relax a structure using Quantum ESPRESSO pw.x."""

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super(PwRelaxWorkChain, cls).define(spec)
        spec.expose_inputs(PwBaseWorkChain, namespace='base',
            exclude=('clean_workdir', 'pw.structure', 'pw.parent_folder'),
            namespace_options={'help': 'Inputs for the `PwBaseWorkChain`.'})
        spec.input('structure', valid_type=orm.StructureData, help='The inputs structure.')
        spec.input('final_scf', valid_type=orm.Bool, default=orm.Bool(False),
            help='If `True`, a final SCF calculation will be performed on the successfully relaxed structure.')
        spec.input('relaxation_scheme', valid_type=orm.Str, default=orm.Str('vc-relax'),
            help='The relaxation scheme to use: choose either `relax` or `vc-relax` for variable cell relax.')
        spec.input('meta_convergence', valid_type=orm.Bool, default=orm.Bool(True),
            help='If `True` the workchain will perform a meta-convergence on the cell volume.')
        spec.input('max_meta_convergence_iterations', valid_type=orm.Int, default=orm.Int(5),
            help='The maximum number of variable cell relax iterations in the meta convergence cycle.')
        spec.input('volume_convergence', valid_type=orm.Float, default=orm.Float(0.01),
            help='The volume difference threshold between two consecutive meta convergence iterations.')
        spec.input('clean_workdir', valid_type=orm.Bool, default=orm.Bool(False),
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
        spec.output('output_structure', valid_type=orm.StructureData, required=True,
            help='The successfully relaxed structure.')

    def setup(self):
        """Input validation and context setup."""
        self.ctx.current_number_of_bands = None
        self.ctx.current_structure = self.inputs.structure
        self.ctx.current_cell_volume = None
        self.ctx.is_converged = False
        self.ctx.iteration = 0

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
        return self.inputs.final_scf.value and self.ctx.is_converged

    def run_relax(self):
        """Run the `PwBaseWorkChain` to run a relax `PwCalculation`."""
        self.ctx.iteration += 1

        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='base'))
        inputs.pw.structure = self.ctx.current_structure
        inputs.pw.parameters = inputs.pw.parameters.get_dict()

        inputs.pw.parameters.setdefault('CONTROL', {})
        inputs.pw.parameters['CONTROL']['calculation'] = self.inputs.relaxation_scheme.value
        inputs.pw.parameters['CONTROL']['restart_mode'] = 'from_scratch'

        # If one of the nested `PwBaseWorkChains` changed the number of bands, apply it here
        if self.ctx.current_number_of_bands is not None:
            inputs.pw.parameters.setdefault('SYSTEM', {})['nbnd'] = self.ctx.current_number_of_bands

        # Set the `CALL` link label
        inputs.metadata.call_link_label = 'iteration_{:02d}'.format(self.ctx.iteration)

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}>'.format(running.pk))

        return ToContext(workchains=append_(running))

    def inspect_relax(self):
        """Inspect the results of the last `PwBaseWorkChain`.

        Compare the cell volume of the relaxed structure of the last completed workchain with the previous. If the
        difference ratio is less than the volume convergence threshold we consider the cell relaxation converged.
        """
        workchain = self.ctx.workchains[-1]

        acceptable_statuses = [
            'ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF'
        ]

        if workchain.is_excepted or workchain.is_killed:
            self.report('relax PwBaseWorkChain was excepted or killed')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        if workchain.is_failed and workchain.exit_status not in PwBaseWorkChain.get_exit_statuses(acceptable_statuses):
            self.report('relax PwBaseWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        try:
            structure = workchain.outputs.output_structure
        except exceptions.NotExistent:
            self.report('relax PwBaseWorkChain finished successful but without output structure')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        prev_cell_volume = self.ctx.current_cell_volume
        curr_cell_volume = structure.get_cell_volume()

        # Set relaxed structure as input structure for next iteration
        self.ctx.current_structure = structure
        self.ctx.current_number_of_bands = workchain.outputs.output_parameters.get_dict()['number_of_bands']
        self.report('after iteration {} cell volume of relaxed structure is {}'
            .format(self.ctx.iteration, curr_cell_volume))

        # After first iteration, simply set the cell volume and restart the next base workchain
        if not prev_cell_volume:
            self.ctx.current_cell_volume = curr_cell_volume

            # If meta convergence is switched off we are done
            if not self.inputs.meta_convergence.value:
                self.ctx.is_converged = True
            return

        # Check whether the cell volume is converged
        volume_threshold = self.inputs.volume_convergence.value
        volume_difference = abs(prev_cell_volume - curr_cell_volume) / prev_cell_volume

        if volume_difference < volume_threshold:
            self.ctx.is_converged = True
            self.report('relative cell volume difference {} smaller than convergence threshold {}'
                .format(volume_difference, volume_threshold))
        else:
            self.report('current relative cell volume difference {} larger than convergence threshold {}'
                .format(volume_difference, volume_threshold))

        self.ctx.current_cell_volume = curr_cell_volume

        return

    def run_final_scf(self):
        """Run the `PwBaseWorkChain` to run a final scf `PwCalculation` for the relaxed structure."""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='base'))
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

        self.report('launching PwBaseWorkChain<{}> for final scf'.format(running.pk))

        return ToContext(workchain_scf=running)

    def inspect_final_scf(self):
        """Inspect the result of the final scf `PwBaseWorkChain`."""
        workchain = self.ctx.workchain_scf

        if not workchain.is_finished_ok:
            self.report('final scf PwBaseWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_FINAL_SCF

    def results(self):
        """Attach the output parameters and structure of the last workchain to the outputs."""
        if self.ctx.is_converged and self.ctx.iteration <= self.inputs.max_meta_convergence_iterations.value:
            self.report('workchain completed after {} iterations'.format(self.ctx.iteration))
        else:
            self.report('maximum number of meta convergence iterations exceeded')

        # Get the latest workchain, which is either the workchain_scf if it ran or otherwise the last regular workchain
        try:
            workchain = self.ctx.workchain_scf
            structure = workchain.inputs.pw__structure
        except AttributeError:
            workchain = self.ctx.workchains[-1]
            structure = workchain.outputs.output_structure

        self.out_many(self.exposed_outputs(workchain, PwBaseWorkChain))
        self.out('output_structure', structure)

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super(PwRelaxWorkChain, self).on_terminated()

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
            self.report('cleaned remote folders of calculations: {}'.format(' '.join(map(str, cleaned_calcs))))
