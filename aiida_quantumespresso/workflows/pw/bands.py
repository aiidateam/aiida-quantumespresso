# -*- coding: utf-8 -*-
from aiida.common.extendeddicts import AttributeDict
from aiida.orm import Code
from aiida.orm.data.base import Str, Float, Bool
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.singlefile import SinglefileData
from aiida.orm.group import Group
from aiida.orm.utils import WorkflowFactory
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, if_
from aiida.work.workfunction import workfunction
from aiida_quantumespresso.utils.mapping import prepare_process_inputs


PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
PwRelaxWorkChain = WorkflowFactory('quantumespresso.pw.relax')


class PwBandsWorkChain(WorkChain):
    """
    Workchain to launch a Quantum Espresso pw.x to calculate a bandstructure for a given
    structure. The structure will first be relaxed followed by a band structure calculation
    """

    @classmethod
    def define(cls, spec):
        super(PwBandsWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        spec.input('pseudo_family', valid_type=Str)
        spec.input('kpoints_distance', valid_type=Float, default=Float(0.2))
        spec.input('kpoints_distance_bands', valid_type=Float, default=Float(0.2))
        spec.input('kpoints', valid_type=KpointsData, required=False)
        spec.input('vdw_table', valid_type=SinglefileData, required=False)
        spec.input('parameters', valid_type=ParameterData)
        spec.input('settings', valid_type=ParameterData, required=False)
        spec.input('options', valid_type=ParameterData, required=False)
        spec.input('automatic_parallelization', valid_type=ParameterData, required=False)
        spec.input('group', valid_type=Str, required=False)
        spec.input_group('relax', required=False)
        spec.outline(
            cls.setup,
            cls.validate_inputs,
            if_(cls.should_do_relax)(
                cls.run_relax,
            ),
            cls.run_seekpath,
            cls.run_scf,
            cls.run_bands,
            cls.results,
        )
        spec.output('primitive_structure', valid_type=StructureData)
        spec.output('seekpath_parameters', valid_type=ParameterData)
        spec.output('scf_parameters', valid_type=ParameterData)
        spec.output('band_parameters', valid_type=ParameterData)
        spec.output('band_structure', valid_type=BandsData)

    def setup(self):
        """
        Initialize context variables that are used during the logical flow of the BaseRestartWorkChain
        """
        self.ctx.inputs = AttributeDict({
            'code': self.inputs.code,
            'pseudo_family': self.inputs.pseudo_family,
            'parameters': self.inputs.parameters.get_dict(),
        })

    def validate_inputs(self):
        """
        Validate inputs that may depend on each other
        """
        if not any([key in self.inputs for key in ['options', 'automatic_parallelization']]):
            self.abort_nowait('you have to specify either the options or automatic_parallelization input')
            return

        # Add the van der Waals kernel table file if specified
        if 'vdw_table' in self.inputs:
            self.ctx.inputs.vdw_table = self.inputs.vdw_table
            self.inputs.relax['vdw_table'] = self.inputs.vdw_table

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

        return

    def should_do_relax(self):
        """
        If the 'relax' input group was specified, we relax the input structure
        """
        return 'relax' in self.inputs

    def run_relax(self):
        """
        Run the PwRelaxWorkChain to run a relax PwCalculation
        """
        inputs = self.inputs.relax
        inputs.update({
            'code': self.inputs.code,
            'structure': self.inputs.structure,
            'pseudo_family': self.inputs.pseudo_family,
        })

        # If options set, add it to the default inputs
        if 'options' in self.inputs:
            inputs['options'] = self.inputs.options

        # If automatic parallelization was set, add it to the default inputs
        if 'automatic_parallelization' in self.inputs:
            inputs['automatic_parallelization'] = self.inputs.automatic_parallelization

        running = submit(PwRelaxWorkChain, **inputs)

        self.report('launching PwRelaxWorkChain<{}>'.format(running.pid))

        return ToContext(workchain_relax=running)

    def run_seekpath(self):
        """
        Run the relaxed structure through SeeKPath to get the new primitive structure, just in case
        the symmetry of the cell changed in the cell relaxation step
        """
        if 'workchain_relax' not in self.ctx:
            structure = self.inputs.structure
        else:
            try:
                structure = self.ctx.workchain_relax.out.output_structure
            except:
                self.abort_nowait('the relax workchain did not output an output_structure node')
                return

        seekpath_parameters = ParameterData(dict={
            'reference_distance': self.inputs.kpoints_distance_bands.value
        })

        result = seekpath_structure_analysis(structure, seekpath_parameters)
        self.ctx.structure = result['primitive_structure']

        # Construct a new kpoint mesh for the scf calculation on the current primitive structure
        kpoints_mesh = KpointsData()
        kpoints_mesh.set_cell_from_structure(self.ctx.structure)
        kpoints_mesh.set_kpoints_mesh_from_density(self.inputs.kpoints_distance.value)

        # Save the kpoints objects for the scf and bands calculation in the context
        self.ctx.kpoints_mesh = kpoints_mesh
        self.ctx.kpoints_path = result['explicit_kpoints']

        self.out('primitive_structure', result['primitive_structure'])
        self.out('seekpath_parameters', result['parameters'])

    def run_scf(self):
        """
        Run the PwBaseWorkChain in scf mode on the primitive cell of (optionally relaxed) input structure
        """
        inputs = self.ctx.inputs
        calculation_mode = 'scf'

        # Set the correct pw.x input parameters
        inputs.parameters['CONTROL']['calculation'] = calculation_mode
        inputs.kpoints = self.ctx.kpoints_mesh
        inputs.structure = self.ctx.structure

        inputs = prepare_process_inputs(inputs)
        running = submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}> in {} mode'.format(running.pid, calculation_mode))

        return ToContext(workchain_scf=running)

    def run_bands(self):
        """
        Run the PwBaseWorkChain to run a bands PwCalculation along the path of high-symmetry determined by Seekpath
        """
        try:
            remote_folder = self.ctx.workchain_scf.out.remote_folder
        except AttributeError as exception:
            self.abort_nowait('the scf workchain did not output a remote_folder node')
            return

        inputs = self.ctx.inputs
        restart_mode = 'restart'
        calculation_mode = 'bands'

        # Set the correct pw.x input parameters
        inputs.parameters['CONTROL']['restart_mode'] = restart_mode
        inputs.parameters['CONTROL']['calculation'] = calculation_mode

        if 'kpoints' in self.inputs:
            inputs.kpoints = self.inputs.kpoints
        else:
            inputs.kpoints = self.ctx.kpoints_path

        inputs.structure = self.ctx.structure
        inputs.parent_folder = remote_folder

        inputs = prepare_process_inputs(inputs)
        running = submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}> in {} mode'.format(running.pid, calculation_mode))

        return ToContext(workchain_bands=running)

    def results(self):
        """
        Attach the desired output nodes directly as outputs of the workchain
        """
        self.report('workchain succesfully completed')
        self.out('scf_parameters', self.ctx.workchain_scf.out.output_parameters)
        self.out('band_parameters', self.ctx.workchain_bands.out.output_parameters)
        self.out('band_structure', self.ctx.workchain_bands.out.output_band)

        if 'group' in self.inputs:
            output_band = self.ctx.workchain_bands.out.output_band
            group, _ = Group.get_or_create(name=self.inputs.group.value)
            group.add_nodes(output_band)
            self.report("storing the output_band<{}> in the group '{}'"
                .format(output_band.pk, self.inputs.group.value))


@workfunction
def seekpath_structure_analysis(structure, parameters):
    """
    This workfunction will take a structure and pass it through SeeKpath to get the
    primitive cell and the path of high symmetry k-points through its Brillouin zone.
    Note that the returned primitive cell may differ from the original structure in
    which case the k-points are only congruent with the primitive cell.
    """
    from aiida.tools import get_explicit_kpoints_path
    return get_explicit_kpoints_path(structure, **parameters.get_dict())
