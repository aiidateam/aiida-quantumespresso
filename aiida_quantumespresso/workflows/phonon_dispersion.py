# -*- coding: utf-8 -*-
"""Workchain to compute the phonon dispersion from the raw initial unrelaxed structure."""
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import ToContext, WorkChain, if_
from aiida.plugins import WorkflowFactory

from aiida_quantumespresso.calculations.functions.seekpath_structure_analysis import seekpath_structure_analysis

PhBaseWorkChain = WorkflowFactory('quantumespresso.ph.base')
PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
Q2rBaseWorkChain = WorkflowFactory('quantumespresso.q2r.base')
MatdynBaseWorkChain = WorkflowFactory('quantumespresso.matdyn.base')


class PhononDispersionWorkChain(WorkChain):
    """Workchain to compute the phonon dispersion from the raw initial unrelaxed structure."""

    defaults = AttributeDict({
        'reference_distance': 0.025,
        'symprec': 1e-05,
        'angle_tolerance': -1.0,
    })

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        super().define(spec)
        spec.expose_inputs(PwBaseWorkChain, namespace='pw')
        spec.expose_inputs(PhBaseWorkChain, namespace='ph', exclude=('parent_folder',))
        spec.expose_inputs(Q2rBaseWorkChain, namespace='q2r', exclude=('parent_folder',))
        spec.expose_inputs(MatdynBaseWorkChain, namespace='matdyn', exclude=('parent_folder',))
        spec.outline(
            cls.setup,
            if_(cls.should_run_pw)(cls.run_pw,),
            if_(cls.should_run_ph)(cls.run_ph,),
            if_(cls.should_run_dispersion)(
                cls.run_q2r,
                cls.run_matdyn,
            ),
            cls.results,
        )
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('output_phonon_bands', valid_type=orm.BandsData)

    def setup(self):
        """Initialize context variables."""
        self.ctx.structure = self.inputs.pw['structure']

    def should_run_pw(self):
        """If a parent_calc is specified in the ph inputs, no PwBaseWorkChain has to be ran."""
        return 'ph' in self.inputs and 'parent_folder' not in self.inputs.ph

    def should_run_ph(self):
        """If a parent_calc is specified in the q2r inputs, no PhBaseWorkChain has to be ran."""
        return 'q2r' in self.inputs and 'parent_folder' not in self.inputs.q2r

    def should_run_dispersion(self):
        """If both the q2r_code and the matdyn_code are passed as inputs, the dispersion is computed."""
        return 'code' in self.inputs.q2r and 'code' in self.inputs.matdyn

    def run_pw(self):
        """Run the PwBaseWorkChain."""
        inputs = AttributeDict(self.inputs.pw)

        running = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launching PwBaseWorkChain<{running.pk}>')

        return ToContext(workflow_pw=running)

    def run_ph(self):
        """Run the PhWorkChain."""
        inputs = AttributeDict(self.inputs.ph)

        # Load parent folder from completed PwBaseWorkChain if not provided in inputs
        if 'parent_folder' not in inputs:
            inputs['parent_folder'] = self.ctx.workflow_pw.out.remote_folder

        running = self.submit(PhBaseWorkChain, **inputs)

        self.report(f'launching PhBaseWorkChain<{running.pk}>')

        return ToContext(workflow_ph=running)

    def run_q2r(self):
        """Run the Q2rCalculation."""
        inputs = AttributeDict(self.inputs.q2r)

        # Load parent folder from completed PhBaseWorkChain if not provided in inputs
        if 'parent_folder' not in inputs:
            inputs['parent_folder'] = self.ctx.workflow_ph.out.retrieved

        running = self.submit(Q2rBaseWorkChain, **inputs)

        self.report(f'launching Q2rBaseWorkChain<{running.pk}>')

        return ToContext(workflow_q2r=running)

    def run_matdyn(self):
        """Run the MatdynCalculation."""
        inputs = AttributeDict(self.inputs.matdyn)

        # Load parent folder from completed PhBaseWorkChain if not provided in inputs
        if 'parent_folder' not in inputs:
            inputs['parent_folder'] = self.ctx.workflow_q2r.out.force_constants

        parameters = {
            'reference_distance': orm.Float(self.defaults.reference_distance),
            'symprec': orm.Float(self.defaults.symprec),
            'angle_tolerance': orm.Float(self.defaults.angle_tolerance),
        }

        # Somehow we need to retrieve the correct structure, which is the one that is used for the
        # pw step. However, if we want to support the skipping of the pw step, we can't obtain the structure from
        # the inputs, but rather will have to obtain it from either the parent_folder input in the ph input group
        # or from the parent folder in the q2r input group if the ph step is also to be skipped. Currently it is
        # not clear how this can be done
        structure = self.ctx.structure
        seekpath_results = seekpath_structure_analysis(structure, **parameters)

        inputs['kpoints'] = seekpath_results['explicit_kpoints']

        running = self.submit(MatdynBaseWorkChain, **inputs)

        self.report(f'launching MatdynBaseWorkChain<{running.pk}>')

        return ToContext(workflow_matdyn=running)

    def results(self):
        """Run the final step after computing the dispersion steps."""
        matdyn_calc = self.ctx.workflow_matdyn

        self.out('output_parameters', matdyn_calc.out.output_parameters)
        self.out('output_phonon_bands', matdyn_calc.out.output_phonon_bands)
