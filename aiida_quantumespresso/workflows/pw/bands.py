# -*- coding: utf-8 -*-
from aiida.orm import Code
from aiida.orm.data.base import Str, Float
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.data.singlefile import SinglefileData
from aiida.orm.utils import WorkflowFactory
from aiida.common.links import LinkType
from aiida.common.exceptions import AiidaException, NotExistent
from aiida.common.datastructures import calc_states
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext, while_, append_
from aiida.work.workfunction import workfunction
from seekpath.aiidawrappers import get_path, get_explicit_k_path

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
        spec.input('kpoints_mesh', valid_type=KpointsData, required=False)
        spec.input('kpoints_distance', valid_type=Float, default=Float(0.2))
        spec.input('vdw_table', valid_type=SinglefileData, required=False)
        spec.input('parameters', valid_type=ParameterData)
        spec.input('settings', valid_type=ParameterData)
        spec.input('options', valid_type=ParameterData)
        spec.input_group('relax')
        spec.outline(
            cls.setup,
            cls.run_relax,
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
        Input validation and context setup
        """
        self.ctx.inputs = {
            'code': self.inputs.code,
            'parameters': self.inputs.parameters.get_dict(),
            'settings': self.inputs.settings,
            'options': self.inputs.options,
        }

        # We expect either a KpointsData with given mesh or a desired distance between k-points
        if all([key not in self.inputs for key in ['kpoints_mesh', 'kpoints_distance']]):
            self.abort_nowait('neither the kpoints_mesh nor a kpoints_distance was specified in the inputs')
            return

        # Add the van der Waals kernel table file if specified
        if 'vdw_table' in self.inputs:
            self.ctx.inputs['vdw_table'] = self.inputs.vdw_table
            self.inputs.relax['vdw_table'] = self.inputs.vdw_table

        # Set the correct relaxation scheme in the input parameters
        if 'CONTROL' not in self.ctx.inputs['parameters']:
            self.ctx.inputs['parameters']['CONTROL'] = {}

        return

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

        running = submit(PwRelaxWorkChain, **inputs)

        self.report('launching PwRelaxWorkChain<{}>'.format(running.pid))

        return ToContext(workchain_relax=running)

    def run_seekpath(self):
        """
        Run the relaxed structure through SeeKPath to get the new primitive structure, just in case
        the symmetry of the cell changed in the cell relaxation step
        """
        try:
            structure = self.ctx.workchain_relax.out.output_structure
        except:
            self.abort_nowait('the relax workchain did not output a output_structure node')
            return

        result = seekpath_structure_analysis(structure)

        self.ctx.structure_relaxed_primitive = result['primitive_structure']
        self.ctx.kpoints_path = result['explicit_kpoints_path']

        self.out('primitive_structure', result['primitive_structure'])
        self.out('seekpath_parameters', result['parameters'])

    def run_scf(self):
        """
        Run the PwBaseWorkChain in scf mode on the primitive cell of the relaxed input structure
        """
        inputs = dict(self.ctx.inputs)
        structure = self.ctx.structure_relaxed_primitive
        calculation_mode = 'scf'

        # Set the correct pw.x input parameters
        inputs['parameters']['CONTROL']['calculation'] = calculation_mode

        # Construct a new kpoint mesh on the current structure or pass the static mesh
        if 'kpoints_distance' in self.inputs:
            kpoints_mesh = KpointsData()
            kpoints_mesh.set_cell_from_structure(structure)
            kpoints_mesh.set_kpoints_mesh_from_density(self.inputs.kpoints_distance.value)
        else:
            kpoints_mesh = self.inputs.kpoints_mesh

        # Final input preparation, wrapping dictionaries in ParameterData nodes
        inputs['kpoints'] = kpoints_mesh
        inputs['structure'] = structure
        inputs['parameters'] = ParameterData(dict=inputs['parameters'])
        inputs['pseudo_family'] = self.inputs.pseudo_family

        running = submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}> in {} mode'.format(running.pid, calculation_mode))

        return ToContext(workchain_scf=running)

    def run_bands(self):
        """
        Run the PwBaseWorkChain to run a bands PwCalculation
        """
        try:
            remote_folder = self.ctx.workchain_scf.out.remote_folder
        except AttributeError as exception:
            self.abort_nowait('the scf workchain did not output a remote_folder node')
            return

        inputs = dict(self.ctx.inputs)
        structure = self.ctx.structure_relaxed_primitive
        restart_mode = 'restart'
        calculation_mode = 'bands'

        # Set the correct pw.x input parameters
        inputs['parameters']['CONTROL']['restart_mode'] = restart_mode
        inputs['parameters']['CONTROL']['calculation'] = calculation_mode

        # Tell the plugin to retrieve the bands
        settings = inputs['settings'].get_dict()
        settings['also_bands'] = True

        # Final input preparation, wrapping dictionaries in ParameterData nodes
        inputs['kpoints'] = self.ctx.kpoints_path
        inputs['structure'] = structure
        inputs['parent_folder'] = remote_folder
        inputs['parameters'] = ParameterData(dict=inputs['parameters'])
        inputs['settings'] = ParameterData(dict=settings)
        inputs['pseudo_family'] = self.inputs.pseudo_family

        running = submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}> in {} mode'.format(running.pid, calculation_mode))

        return ToContext(workchain_bands=running)

    def results(self):
        """
        Attach the desired output nodes directly as outputs of the workchain
        """
        calculation_band = self.ctx.workchain_bands.get_outputs(link_type=LinkType.CALL)[0]

        self.report('workchain succesfully completed')
        self.out('scf_parameters', self.ctx.workchain_scf.out.output_parameters)
        self.out('band_parameters', calculation_band.out.output_parameters)
        self.out('band_structure', calculation_band.out.output_band)


@workfunction
def seekpath_structure_analysis(structure):
    """
    This workfunction will take a structure and pass it through SeeKpath to get the
    primitive cell and the path of high symmetry k-points through its Brillouin zone.
    Note that the returned primitive cell may differ from the original structure in
    which case the k-points are only congruent with the primitive cell.
    """
    seekpath_info = get_path(structure)
    explicit_path = get_explicit_k_path(structure)

    primitive_structure = seekpath_info.pop('primitive_structure')
    conv_structure = seekpath_info.pop('conv_structure')
    parameters = ParameterData(dict=seekpath_info)

    result = {
        'parameters': parameters,
        'conv_structure': conv_structure,
        'primitive_structure': primitive_structure,
        'explicit_kpoints_path': explicit_path['explicit_kpoints'],
    }

    return result