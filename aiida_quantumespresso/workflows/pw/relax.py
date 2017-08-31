# -*- coding: utf-8 -*-
from aiida.orm import Code
from aiida.orm.data.base import Str, Float
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.singlefile import SinglefileData
from aiida.common.exceptions import AiidaException, NotExistent
from aiida.common.datastructures import calc_states
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, while_, append_
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

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
        spec.input('vdw_table', valid_type=SinglefileData, required=False)
        spec.input('parameters', valid_type=ParameterData)
        spec.input('settings', valid_type=ParameterData)
        spec.input('options', valid_type=ParameterData)
        spec.input('relaxation_scheme', valid_type=Str, default=Str('vc-relax'))
        spec.input('volume_convergence', valid_type=Float, default=Float(0.01))
        spec.outline(
            cls.setup,
            while_(cls.should_run_relax)(
                cls.run_relax,
                cls.inspect_relax,
            ),
            cls.results,
        )
        spec.dynamic_output()

    def setup(self):
        """
        Input validation and context setup
        """
        self.ctx.current_cell_volume = None
        self.ctx.is_converged = False
        self.ctx.iteration = 0

        self.ctx.inputs = {
            'code': self.inputs.code,
            'structure': self.inputs.structure,
            'parameters': self.inputs.parameters.get_dict(),
            'settings': self.inputs.settings,
            'options': self.inputs.options,
        }

        # We expect either a KpointsData with given mesh or a desired distance between k-points
        if all([key not in self.inputs for key in ['kpoints', 'kpoints_distance']]):
            self.abort_nowait('neither the kpoints nor a kpoints_distance was specified in the inputs')
            return

        # We expect either a pseudo family string or an explicit list of pseudos
        if self.inputs.pseudo_family:
            self.ctx.inputs['pseudo_family'] = self.inputs.pseudo_family
        elif self.inputs.pseudos:
            self.ctx.inputs['pseudos'] = self.inputs.pseudos
        else:
            self.abort_nowait('neither explicit pseudos nor a pseudo_family was specified in the inputs')
            return

        # Add the van der Waals kernel table file if specified
        if 'vdw_table' in self.inputs:
            self.ctx.inputs['vdw_table'] = self.inputs.vdw_table

        # Set the correct relaxation scheme in the input parameters
        if 'CONTROL' not in self.ctx.inputs['parameters']:
            self.ctx.inputs['parameters']['CONTROL'] = {}

        self.ctx.inputs['parameters']['CONTROL']['calculation'] = self.inputs.relaxation_scheme

        return

    def should_run_relax(self):
        """
        Return whether a relaxation workchain should be run, which is the case as long as the volume
        change between two consecutive relaxation runs is larger than the specified volumen convergence
        threshold value.
        """
        return not self.ctx.is_converged

    def run_relax(self):
        """
        Run the PwBaseWorkChain to run a relax PwCalculation
        """
        self.ctx.iteration += 1

        inputs = dict(self.ctx.inputs)
        inputs['parameters'] = ParameterData(dict=inputs['parameters'])

        # Construct a new kpoint mesh on the current structure or pass the static mesh
        if 'kpoints' not in self.inputs or self.inputs.kpoints == None:
            kpoints = KpointsData()
            kpoints.set_cell_from_structure(self.ctx.inputs['structure'])
            kpoints.set_kpoints_mesh_from_density(self.inputs.kpoints_distance.value, force_parity=True)
            inputs['kpoints'] = kpoints
        else:
            inputs['kpoints'] = self.inputs.kpoints

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
        self.ctx.inputs['structure'] = structure
        self.report('after iteration {} cell volume of relaxed structure is {}'
            .format(self.ctx.iteration, curr_cell_volume))

        # After first iteration, simply set the cell volume and restart the next base workchain
        if not prev_cell_volume:
            self.ctx.current_cell_volume = curr_cell_volume
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

    def results(self):
        """
        Attach the output parameters and structure of the last workchain to the outputs
        """
        self.report('workchain completed after {} iterations'.format(self.ctx.iteration))

        workchain = self.ctx.workchains[-1]
        link_labels = ['output_structure', 'output_parameters', 'remote_folder']

        for link_label in link_labels:
            if link_label in workchain.out:
                node = workchain.get_outputs_dict()[link_label]
                self.out(link_label, node)
                self.report("attaching {}<{}> as an output node with label '{}'"
                    .format(node.__class__.__name__, node.pk, link_label))