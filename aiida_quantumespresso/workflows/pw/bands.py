# -*- coding: utf-8 -*-
from aiida.common.extendeddicts import AttributeDict
from aiida.orm.calculation import JobCalculation
from aiida.orm.data.base import Str, Bool
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.group import Group
from aiida.orm.utils import WorkflowFactory
from aiida.work.workchain import WorkChain, ToContext, if_
from aiida.work.workfunctions import workfunction
from aiida_quantumespresso.utils.mapping import prepare_process_inputs


PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
PwRelaxWorkChain = WorkflowFactory('quantumespresso.pw.relax')


class PwBandsWorkChain(WorkChain):
    """Workchain to compute a band structure for a given structure using Quantum ESPRESSO pw.x"""

    @classmethod
    def define(cls, spec):
        super(PwBandsWorkChain, cls).define(spec)
        spec.expose_inputs(PwRelaxWorkChain, namespace='relax', exclude=('structure',))
        spec.expose_inputs(PwBaseWorkChain, namespace='scf', exclude=('structure', 'kpoints'))
        spec.expose_inputs(PwBaseWorkChain, namespace='bands', exclude=('structure',))
        spec.input('structure', valid_type=StructureData)
        spec.input('clean_workdir', valid_type=Bool, default=Bool(False))
        spec.input('group', valid_type=Str, required=False)
        spec.outline(
            cls.setup,
            if_(cls.should_do_relax)(
                cls.run_relax,
                cls.inspect_relax,
            ),
            cls.run_seekpath,
            cls.run_scf,
            cls.inspect_scf,
            cls.run_bands,
            cls.inspect_bands,
            cls.results,
        )
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_RELAX',
            message='the PwRelaxWorkChain sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_SCF',
            message='the scf PwBasexWorkChain sub process failed')
        spec.exit_code(403, 'ERROR_SUB_PROCESS_FAILED_BANDS',
            message='the bands PwBasexWorkChain sub process failed')
        spec.output('primitive_structure', valid_type=StructureData)
        spec.output('seekpath_parameters', valid_type=ParameterData)
        spec.output('scf_parameters', valid_type=ParameterData)
        spec.output('band_parameters', valid_type=ParameterData)
        spec.output('band_structure', valid_type=BandsData)

    def setup(self):
        """Define the current structure in the context to be the input structure."""
        self.ctx.current_structure = self.inputs.structure

    def should_do_relax(self):
        """If the 'relax' input namespace was specified, we relax the input structure."""
        return 'relax' in self.inputs

    def run_relax(self):
        """Run the PwRelaxWorkChain to run a relax PwCalculation."""
        inputs = AttributeDict(self.exposed_inputs(PwRelaxWorkChain, namespace='relax'))
        inputs.structure = self.ctx.current_structure

        running = self.submit(PwRelaxWorkChain, **inputs)

        self.report('launching PwRelaxWorkChain<{}>'.format(running.pk))

        return ToContext(workchain_relax=running)

    def inspect_relax(self):
        """Verify that the PwRelaxWorkChain finished successfully."""
        workchain = self.ctx.workchain_relax

        if not workchain.is_finished_ok:
            self.report('PwRelaxWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX
        else:
            self.ctx.current_structure = workchain.out.output_structure

    def run_seekpath(self):
        """
        Run the relaxed structure through SeeKPath to get the new primitive structure, just in case
        the symmetry of the cell changed in the cell relaxation step
        """
        if 'kpoints_distance' in self.inputs.bands:
            seekpath_parameters = ParameterData(dict={
                'reference_distance': self.inputs.bands.kpoints_distance.value
            })
        else:
            seekpath_parameters = ParameterData(dict={})

        result = seekpath_structure_analysis(self.ctx.current_structure, seekpath_parameters)
        self.ctx.current_structure = result['primitive_structure']
        self.ctx.kpoints_path = result['explicit_kpoints']

        self.out('primitive_structure', result['primitive_structure'])
        self.out('seekpath_parameters', result['parameters'])

    def run_scf(self):
        """Run the PwBaseWorkChain in scf mode on the primitive cell of (optionally relaxed) input structure"""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        inputs.structure = self.ctx.current_structure
        inputs.parameters = inputs.parameters.get_dict()
        inputs.parameters.setdefault('CONTROL', {})
        inputs.parameters['CONTROL']['calculation'] = 'scf'

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}> in {} mode'.format(running.pk, 'scf'))

        return ToContext(workchain_scf=running)

    def inspect_scf(self):
        """Verify that the PwBaseWorkChain for the scf run finished successfully."""
        workchain = self.ctx.workchain_bands

        if not workchain.is_finished_ok:
            self.report('scf PwBaseWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF
        else:
            self.ctx.current_folder = workchain.out.remote_folder

    def run_bands(self):
        """Run the PwBaseWorkChain in bands mode along the path of high-symmetry determined by Seekpath."""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='bands'))
        inputs.parameters = inputs.parameters.get_dict()
        inputs.parameters.setdefault('CONTROL', {})
        inputs.parameters['CONTROL']['restart_mode'] = 'restart'
        inputs.parameters['CONTROL']['calculation'] = 'bands'

        if 'kpoints' not in self.inputs.bands:
            inputs.kpoints = self.ctx.kpoints_path

        inputs.structure = self.ctx.current_structure
        inputs.parent_folder = self.ctx.current_folder

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}> in {} mode'.format(running.pk, 'bands'))

        return ToContext(workchain_bands=running)

    def inspect_bands(self):
        """Verify that the PwBaseWorkChain for the bands run finished successfully."""
        workchain = self.ctx.workchain_bands

        if not workchain.is_finished_ok:
            self.report('bands PwBaseWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_BANDS

    def results(self):
        """Attach the desired output nodes directly as outputs of the workchain."""
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

    def on_terminated(self):
        """
        If the clean_workdir input was set to True, recursively collect all called Calculations by
        ourselves and our called descendants, and clean the remote folder for the JobCalculation instances
        """
        super(PwBandsWorkChain, self).on_terminated()

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
