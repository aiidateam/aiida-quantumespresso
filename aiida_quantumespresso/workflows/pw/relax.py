# -*- coding: utf-8 -*-
from aiida.common.links import LinkType
from aiida.orm import Code
from aiida.orm.data.base import Bool, Float, Int, Str 
from aiida.orm.data.folder import FolderData
from aiida.orm.data.remote import RemoteData
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.singlefile import SinglefileData
from aiida.orm.group import Group
from aiida.orm.utils import CalculationFactory, WorkflowFactory
from aiida.common.extendeddicts import AttributeDict
from aiida.common.exceptions import AiidaException, NotExistent
from aiida.common.datastructures import calc_states
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, if_, while_, append_
from aiida_quantumespresso.utils.mapping import prepare_process_inputs


PwCalculation = CalculationFactory('quantumespresso.pw')
PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')


class PwRelaxWorkChain(WorkChain):
    """
    Workchain to launch a Quantum Espresso pw.x to relax a structure
    """

    @classmethod
    def define(cls, spec):
        super(PwRelaxWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        spec.input_group('pseudos', required=False)
        spec.input('pseudo_family', valid_type=Str, required=False)
        spec.input('kpoints', valid_type=KpointsData, required=False)
        spec.input('kpoints_distance', valid_type=Float, default=Float(0.2))
        spec.input('kpoints_force_parity', valid_type=Bool, default=Bool(False))
        spec.input('vdw_table', valid_type=SinglefileData, required=False)
        spec.input('parameters', valid_type=ParameterData)
        spec.input('settings', valid_type=ParameterData, required=False)
        spec.input('options', valid_type=ParameterData, required=False)
        spec.input('automatic_parallelization', valid_type=ParameterData, required=False)
        spec.input('final_scf', valid_type=Bool, default=Bool(False))
        spec.input('group', valid_type=Str, required=False)
        spec.input('max_iterations', valid_type=Int, default=Int(5))
        spec.input('max_meta_convergence_iterations', valid_type=Int, default=Int(5))
        spec.input('meta_convergence', valid_type=Bool, default=Bool(True))
        spec.input('relaxation_scheme', valid_type=Str, default=Str('vc-relax'))
        spec.input('volume_convergence', valid_type=Float, default=Float(0.01))
        spec.input('clean_workdir', valid_type=Bool, default=Bool(False))
        spec.outline(
            cls.setup,
            cls.validate_inputs,
            while_(cls.should_run_relax)(
                cls.run_relax,
                cls.inspect_relax,
            ),
            if_(cls.should_run_final_scf)(
                cls.run_final_scf,
            ),
            cls.results,
        )
        spec.output('output_structure', valid_type=StructureData)
        spec.output('output_parameters', valid_type=ParameterData)
        spec.output('remote_folder', valid_type=RemoteData)
        spec.output('retrieved', valid_type=FolderData)

    def setup(self):
        """
        Input validation and context setup
        """
        self.ctx.current_structure = self.inputs.structure
        self.ctx.current_parent_folder = None
        self.ctx.current_cell_volume = None
        self.ctx.is_converged = False
        self.ctx.iteration = 0

        self.ctx.inputs = AttributeDict({
            'code': self.inputs.code,
            'parameters': self.inputs.parameters.get_dict(),
            'max_iterations': self.inputs.max_iterations
        })

        # We expect either a KpointsData with given mesh or a desired distance between k-points
        if all([key not in self.inputs for key in ['kpoints', 'kpoints_distance']]):
            self.abort_nowait('neither the kpoints nor a kpoints_distance was specified in the inputs')
            return

        # We expect either a pseudo family string or an explicit list of pseudos
        if self.inputs.pseudo_family:
            self.ctx.inputs.pseudo_family = self.inputs.pseudo_family
        elif self.inputs.pseudos:
            self.ctx.inputs.pseudos = self.inputs.pseudos
        else:
            self.abort_nowait('neither explicit pseudos nor a pseudo_family was specified in the inputs')
            return

        # Add the van der Waals kernel table file if specified
        if 'vdw_table' in self.inputs:
            self.ctx.inputs.vdw_table = self.inputs.vdw_table

        # Set the correct relaxation scheme in the input parameters
        if 'CONTROL' not in self.ctx.inputs.parameters:
            self.ctx.inputs.parameters['CONTROL'] = {}

        if 'settings' in self.inputs:
            self.ctx.inputs.settings = self.inputs.settings

        # If options set, add it to the default inputs
        if 'options' in self.inputs:
            self.ctx.inputs.options = self.inputs.options

        # If automatic parallelization was set, add it to the default inputs
        if 'automatic_parallelization' in self.inputs:
            self.ctx.inputs.automatic_parallelization = self.inputs.automatic_parallelization

        self.ctx.inputs.parameters['CONTROL']['calculation'] = self.inputs.relaxation_scheme.value

        return

    def validate_inputs(self):
        """
        Validate inputs that may depend on each other
        """
        if not any([key in self.inputs for key in ['options', 'automatic_parallelization']]):
            self.abort_nowait('you have to specify either the options or automatic_parallelization input')
            return

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

        inputs = self.ctx.inputs
        inputs['structure'] = self.ctx.current_structure

        # Construct a new kpoint mesh on the current structure or pass the static mesh
        if 'kpoints' not in self.inputs or self.inputs.kpoints == None:
            kpoints = KpointsData()
            kpoints.set_cell_from_structure(self.ctx.current_structure)
            kpoints.set_kpoints_mesh_from_density(
                self.inputs.kpoints_distance.value,
                force_parity=self.inputs.kpoints_force_parity.value
            )
            inputs['kpoints'] = kpoints
        else:
            inputs['kpoints'] = self.inputs.kpoints

        inputs = prepare_process_inputs(inputs)
        running = submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}>'.format(running.pid))

        return ToContext(workchains=append_(running))

    def inspect_relax(self):
        """
        Compare the cell volume of the relaxed structure of the last completed workchain with the previous.
        If the difference ratio is less than the volume convergence threshold we consider the cell relaxation
        converged and can quit the workchain. If the 
        """
        try:
            workchain = self.ctx.workchains[-1]
        except IndexError:
            self.abort_nowait('the first iteration finished without returning a PwBaseWorkChain')
            return

        try:
            structure = workchain.out.output_structure
        except AttributeError as exception:
            self.abort_nowait('the workchain did not have an output structure and probably failed')
            return

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
        """
        Run the PwBaseWorkChain to run a final scf PwCalculation for the relaxed structure
        """
        inputs = self.ctx.inputs
        structure = self.ctx.current_structure

        inputs.parameters['CONTROL']['calculation'] = 'scf'
        inputs.parameters['CONTROL']['restart_mode'] = 'restart'

        # Construct a new kpoint mesh on the current structure or pass the static mesh
        if 'kpoints' not in self.inputs or self.inputs.kpoints == None:
            kpoints = KpointsData()
            kpoints.set_cell_from_structure(structure)
            kpoints.set_kpoints_mesh_from_density(
                self.inputs.kpoints_distance.value,
                force_parity=self.inputs.kpoints_force_parity.value
            )
        else:
            kpoints = self.inputs.kpoints

        inputs.update({
            'kpoints': kpoints,
            'structure': structure,
            'parent_folder': self.ctx.current_parent_folder,
        })

        inputs = prepare_process_inputs(inputs)
        running = submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}> for final scf'.format(running.pid))

        return ToContext(workchain_scf=running)

    def results(self):
        """
        Attach the output parameters and structure of the last workchain to the outputs
        """
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
            except (AttributeError, IndexError) as exception:
                self.report('could not retrieve the last run PwCalculation to add to the result group')
            else:
                group, _ = Group.get_or_create(name=self.inputs.group.value)
                group.add_nodes(calculation)
                self.report("storing the final PwCalculation<{}> in the group '{}'"
                    .format(calculation.pk, self.inputs.group.value))

        # Store the structure separately since the PwBaseWorkChain of final scf will not have it as an output
        self.out('output_structure', structure)
        self.report("attaching {}<{}> as an output node with label '{}'"
            .format(structure.__class__.__name__, structure.pk, 'output_structure'))

        for link_label in ['output_parameters', 'remote_folder',  'retrieved']:
            if link_label in workchain.out:
                node = workchain.get_outputs_dict()[link_label]
                self.out(link_label, node)
                self.report("attaching {}<{}> as an output node with label '{}'"
                    .format(node.__class__.__name__, node.pk, link_label))

    def on_destroy(self):
        """
        Clean remote folders of all PwCalculations run by ourselves and the called subworkchains, if the clean_workdir
        parameter was set to true in the Workchain inputs. We perform this cleaning only at the very end of the
        workchain and do not pass the clean_workdir input directly to the sub workchains that we call, because some of
        the sub workchains may rely on the calculation of on of the previous sub workchains.
        """
        super(PwRelaxWorkChain, self).on_destroy()
        if not self.has_finished():
            return

        if not self.inputs.clean_workdir.value:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = []

        try:
            workchains = self.ctx.workchains
        except AttributeError:
            workchains = []

        try:
            workchains.append(self.ctx.workchain_scf)
        except AttributeError:
            pass

        for workchain in workchains:
            calculations = workchain.get_outputs(link_type=LinkType.CALL)

            for calculation in calculations:
                if isinstance(calculation, PwCalculation):
                    try:
                        calculation.out.remote_folder._clean()
                        cleaned_calcs.append(calculation.pk)
                    except Exception:
                        pass

        if len(cleaned_calcs) > 0:
            self.report('cleaned remote folders of calculations: {}'.format(' '.join(map(str, cleaned_calcs))))
