# -*- coding: utf-8 -*-
"""Workchain to compute a band structure for a given structure using Quantum ESPRESSO pw.x"""
from __future__ import absolute_import

from six.moves import map

from aiida import orm
from aiida.common import AttributeDict
from aiida.plugins import WorkflowFactory
from aiida.engine import WorkChain, ToContext, if_

from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.workflows.functions.seekpath_structure_analysis import seekpath_structure_analysis

PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
PwRelaxWorkChain = WorkflowFactory('quantumespresso.pw.relax')


class PwBandsWorkChain(WorkChain):
    """Workchain to compute a band structure for a given structure using Quantum ESPRESSO pw.x"""

    @classmethod
    def define(cls, spec):
        # yapf: disable
        super(PwBandsWorkChain, cls).define(spec)
        spec.expose_inputs(PwRelaxWorkChain, namespace='relax', exclude=('clean_workdir', 'structure'),
            namespace_options={'required': False, 'populate_defaults': False,
            'help': 'Inputs for the `PwRelaxWorkChain`, if not specified at all, the relaxation step is skipped.'})
        spec.expose_inputs(PwBaseWorkChain, namespace='scf', exclude=('clean_workdir', 'pw.structure'),
            namespace_options={'help': 'Inputs for the `PwBaseWorkChain` for the SCF calculation.'})
        spec.expose_inputs(PwBaseWorkChain, namespace='bands', exclude=('clean_workdir', 'pw.structure'),
            namespace_options={'help': 'Inputs for the `PwBaseWorkChain` for the BANDS calculation.'})
        spec.input('structure', valid_type=orm.StructureData, help='The inputs structure.')
        spec.input('clean_workdir', valid_type=orm.Bool, default=orm.Bool(False),
            help='If `True`, work directories of all called calculation will be cleaned at the end of execution.')
        spec.input('nbands_factor', valid_type=orm.Float, default=orm.Float(1.2),
            help='The number of bands for the BANDS calculation is that used for the SCF multiplied by this factor.')
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
        spec.output('primitive_structure', valid_type=orm.StructureData,
            help='The normalized and primitivized structure for which the bands are computed.')
        spec.output('seekpath_parameters', valid_type=orm.Dict,
            help='The parameters used in the SeeKpath call to normalize the input or relaxed structure.')
        spec.output('scf_parameters', valid_type=orm.Dict,
            help='The output parameters of the SCF `PwBaseWorkChain`.')
        spec.output('band_parameters', valid_type=orm.Dict,
            help='The output parameters of the BANDS `PwBaseWorkChain`.')
        spec.output('band_structure', valid_type=orm.BandsData,
            help='The computed band structure.')

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

        self.ctx.current_structure = workchain.outputs.output_structure

    def run_seekpath(self):
        """Run the structure through SeeKpath to get the primitive and normalized structure.

        This is performed regardless of whether the inputs structure was relaxed.
        """
        if 'kpoints_distance' in self.inputs.bands:
            seekpath_parameters = orm.Dict(dict={
                'reference_distance': self.inputs.bands.kpoints_distance.value
            })
        else:
            seekpath_parameters = orm.Dict(dict={})

        result = seekpath_structure_analysis(self.ctx.current_structure, seekpath_parameters)
        self.ctx.current_structure = result['primitive_structure']
        self.ctx.kpoints_path = result['explicit_kpoints']

        self.out('primitive_structure', result['primitive_structure'])
        self.out('seekpath_parameters', result['parameters'])

    def run_scf(self):
        """Run the PwBaseWorkChain in scf mode on the primitive cell of (optionally relaxed) input structure"""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        inputs.pw.structure = self.ctx.current_structure
        inputs.pw.parameters = inputs.pw.parameters.get_dict()
        inputs.pw.parameters.setdefault('CONTROL', {})
        inputs.pw.parameters['CONTROL']['calculation'] = 'scf'

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}> in {} mode'.format(running.pk, 'scf'))

        return ToContext(workchain_scf=running)

    def inspect_scf(self):
        """Verify that the PwBaseWorkChain for the scf run finished successfully."""
        workchain = self.ctx.workchain_scf

        if not workchain.is_finished_ok:
            self.report('scf PwBaseWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF

        self.ctx.current_folder = workchain.outputs.remote_folder

    def run_bands(self):
        """Run the PwBaseWorkChain in bands mode along the path of high-symmetry determined by seekpath."""

        # Get info from SCF on number of electrons and number of spin components
        scf_out_dict = self.ctx.workchain_scf.outputs.output_parameters.get_dict()
        nelectron = int(scf_out_dict['number_of_electrons'])
        nspin = int(scf_out_dict['number_of_spin_components'])
        nbands = max(
            int(0.5 * nelectron * nspin * self.inputs.nbands_factor.value),
            int(0.5 * nelectron * nspin) + 4 * nspin)

        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='bands'))
        inputs.pw.parameters = inputs.pw.parameters.get_dict()

        inputs.pw.parameters.setdefault('CONTROL', {})
        inputs.pw.parameters.setdefault('SYSTEM', {})
        inputs.pw.parameters.setdefault('ELECTRONS', {})

        inputs.pw.parameters['CONTROL']['restart_mode'] = 'restart'
        inputs.pw.parameters['CONTROL']['calculation'] = 'bands'
        inputs.pw.parameters['ELECTRONS']['diagonalization'] = 'cg'
        inputs.pw.parameters['ELECTRONS']['diago_full_acc'] = True
        inputs.pw.parameters['SYSTEM']['nbnd'] = nbands

        if 'kpoints' not in self.inputs.bands:
            inputs.kpoints = self.ctx.kpoints_path

        inputs.pw.structure = self.ctx.current_structure
        inputs.pw.parent_folder = self.ctx.current_folder

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
        self.out('scf_parameters', self.ctx.workchain_scf.outputs.output_parameters)
        self.out('band_parameters', self.ctx.workchain_bands.outputs.output_parameters)
        self.out('band_structure', self.ctx.workchain_bands.outputs.output_band)

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super(PwBandsWorkChain, self).on_terminated()

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
