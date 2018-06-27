# -*- coding: utf-8 -*-
from aiida.common.extendeddicts import AttributeDict
from aiida.orm.calculation import JobCalculation
from aiida.orm.data.base import Bool, Float, Int, Str
from aiida.orm.data.structure import StructureData
from aiida.orm.group import Group
from aiida.orm.utils import CalculationFactory, WorkflowFactory
from aiida.work.workchain import WorkChain, ToContext, if_, while_, append_
from aiida_quantumespresso.utils.mapping import prepare_process_inputs


PwCalculation = CalculationFactory('quantumespresso.pw')
PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')


class PwRelaxWorkChain(WorkChain):
    """Workchain to relax a structure using Quantum ESPRESSO pw.x"""

    @classmethod
    def define(cls, spec):
        super(PwRelaxWorkChain, cls).define(spec)
        spec.expose_inputs(PwBaseWorkChain, namespace='base', exclude=('structure',))
        spec.input('structure', valid_type=StructureData)
        spec.input('final_scf', valid_type=Bool, default=Bool(False))
        spec.input('group', valid_type=Str, required=False)
        spec.input('relaxation_scheme', valid_type=Str, default=Str('vc-relax'))
        spec.input('meta_convergence', valid_type=Bool, default=Bool(True))
        spec.input('max_meta_convergence_iterations', valid_type=Int, default=Int(5))
        spec.input('volume_convergence', valid_type=Float, default=Float(0.01))
        spec.input('clean_workdir', valid_type=Bool, default=Bool(False))
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
        spec.expose_outputs(PwBaseWorkChain)

    def setup(self):
        """
        Input validation and context setup
        """
        self.ctx.current_structure = self.inputs.structure
        self.ctx.current_parent_folder = None
        self.ctx.current_cell_volume = None
        self.ctx.is_converged = False
        self.ctx.iteration = 0

    def should_run_relax(self):
        """
        Return whether a relaxation workchain should be run, which is the case as long as the volume
        change between two consecutive relaxation runs is larger than the specified volume convergence
        threshold value and the maximum number of meta convergence iterations is not exceeded
        """
        return not self.ctx.is_converged and self.ctx.iteration < self.inputs.max_meta_convergence_iterations.value

    def should_run_final_scf(self):
        """
        Return whether after successful relaxation a final scf calculation should be run. If the maximum number of
        meta convergence iterations has been exceeded and convergence has not been reached, the structure cannot be
        considered to be relaxed and the final scf should not be run
        """
        return self.inputs.final_scf.value and self.ctx.is_converged

    def run_relax(self):
        """
        Run the PwBaseWorkChain to run a relax PwCalculation
        """
        self.ctx.iteration += 1

        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='base'))
        inputs.structure = self.ctx.current_structure
        inputs.parameters = inputs.parameters.get_dict()

        inputs.parameters.setdefault('CONTROL', {})
        inputs.parameters['CONTROL']['calculation'] = self.inputs.relaxation_scheme.value

        # Do not clean workdirs of sub workchains, because then we won't be able to restart from them
        inputs.pop('clean_workdir', None)

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}>'.format(running.pk))

        return ToContext(workchains=append_(running))

    def inspect_relax(self):
        """
        Compare the cell volume of the relaxed structure of the last completed workchain with the previous.
        If the difference ratio is less than the volume convergence threshold we consider the cell relaxation
        converged and can quit the workchain. If the
        """
        workchain = self.ctx.workchains[-1]

        if not self.workchain.is_finished_ok:
            self.report('relax PwBaseWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX
        else:
            structure = workchain.out.output_structure

        prev_cell_volume = self.ctx.current_cell_volume
        curr_cell_volume = structure.get_cell_volume()

        # Set relaxed structure as input structure for next iteration
        self.ctx.current_parent_folder = workchain.out.remote_folder
        self.ctx.current_structure = structure
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
        """Run the PwBaseWorkChain to run a final scf PwCalculation for the relaxed structure."""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='base'))
        inputs.structure = self.ctx.current_structure
        inputs.parent_folder = self.ctx.current_parent_folder
        inputs.parameters = inputs.parameters.get_dict()

        inputs.parameters.setdefault('CONTROL', {})
        inputs.parameters['CONTROL']['calculation'] = 'scf'
        inputs.parameters['CONTROL']['restart_mode'] = 'restart'

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}> for final scf'.format(running.pk))

        return ToContext(workchain_scf=running)

    def inspect_final_scf(self):
        """Inspect the result of the final scf PwBaseWorkChain."""
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
            structure = workchain.inp.structure
        except AttributeError:
            workchain = self.ctx.workchains[-1]
            structure = workchain.out.output_structure

        if 'group' in self.inputs:
            # Retrieve the final successful PwCalculation through the output_parameters of the PwBaseWorkChain
            try:
                calculation = workchain.out.output_parameters.get_inputs(node_type=PwCalculation)[0]
            except (AttributeError, IndexError):
                self.report('could not retrieve the last run PwCalculation to add to the result group')
            else:
                group, _ = Group.get_or_create(name=self.inputs.group.value)
                group.add_nodes(calculation)
                self.report("storing the final PwCalculation<{}> in the group '{}'"
                    .format(calculation.pk, self.inputs.group.value))

        self.out_many(self.exposed_outputs(workchain, PwBaseWorkChain))
        self.out('output_structure', structure)

    def on_terminated(self):
        """
        If the clean_workdir input was set to True, recursively collect all called Calculations by
        ourselves and our called descendants, and clean the remote folder for the JobCalculation instances
        """
        super(PwRelaxWorkChain, self).on_terminated()

        if self.inputs.clean_workdir.value is False:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = []

        for called_descendant in self.calc.called_descendants:
            if isinstance(called_descendant, JobCalculation):
                try:
                    called_descendant.out.remote_folder._clean()
                    cleaned_calcs.append(called_descendant.pk)
                except (IOError, OSError, KeyError):
                    pass

        if cleaned_calcs:
            self.report('cleaned remote folders of calculations: {}'.format(' '.join(map(str, cleaned_calcs))))
