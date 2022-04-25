# -*- coding: utf-8 -*-
"""Workchain to perform a ph.x calculation with optional parallelization over q-points."""
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import ToContext, WorkChain, if_
from aiida.plugins import CalculationFactory, WorkflowFactory

PwCalculation = CalculationFactory('quantumespresso.pw')
PhBaseWorkChain = WorkflowFactory('quantumespresso.ph.base')
PhParallelizeQpointsWorkChain = WorkflowFactory('quantumespresso.ph.parallelize_qpoints')


class PhWorkChain(WorkChain):
    """Workchain that will run a Quantum Espresso ph.x calculation based on a previously completed pw.x calculation.

    If specified through the 'parallelize_qpoints' boolean input parameter, the calculation will be parallelized over
    the provided q-points by running the `PhParallelizeQpointsWorkChain`. Otherwise a single `PhBaseWorkChain` will be
    launched that will compute every q-point serially.
    """

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        super().define(spec)
        spec.input('code', valid_type=orm.Code)
        spec.input('parent_folder', valid_type=orm.RemoteData)
        spec.input('qpoints', valid_type=orm.KpointsData)
        spec.input('parameters', valid_type=orm.Dict, required=False)
        spec.input('settings', valid_type=orm.Dict, required=False)
        spec.input('options', valid_type=orm.Dict, required=False)
        spec.input('clean_workdir', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.input('compute_epsil', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.input('alpha_mix', valid_type=orm.Float, default=lambda: orm.Float(0.7))
        spec.input('max_iterations', valid_type=orm.Int, default=lambda: orm.Int(10))
        spec.input('parallelize_qpoints', valid_type=orm.Bool, default=lambda: orm.Bool(True))
        spec.outline(
            cls.setup,
            if_(cls.should_run_parallel)(cls.run_parallel,).else_(
                cls.run_serial,
            ),
            cls.results,
        )
        spec.output('retrieved', valid_type=orm.FolderData)

    def setup(self):
        """Initialize and validate the inputs."""
        self.ctx.inputs = AttributeDict({
            'code': self.inputs.code,
            'parent_folder': self.inputs.parent_folder,
            'qpoints': self.inputs.qpoints,
            'clean_workdir': self.inputs.clean_workdir,
            'compute_epsil': self.inputs.compute_epsil,
            'alpha_mix': self.inputs.alpha_mix,
            'max_iterations': self.inputs.max_iterations,
        })

        if 'parameters' in self.inputs:
            self.ctx.inputs['parameters'] = self.inputs.parameters

        if 'settings' in self.inputs:
            self.ctx.inputs['settings'] = self.inputs.settings

        if 'options' in self.inputs:
            self.ctx.inputs['options'] = self.inputs.options

    def should_run_parallel(self):
        """Return whether the calculation should be parallelized over the qpoints."""
        return self.inputs.parallelize_qpoints

    def run_parallel(self):
        """Run the PhParallelizeQpointsWorkChain."""
        running = self.submit(PhParallelizeQpointsWorkChain, **self.ctx.inputs)
        self.report(f'running in parallel, launching PhParallelizeQpointsWorkChain<{running.pk}>')
        return ToContext(workchain=running)

    def run_serial(self):
        """Run the PhBaseWorkChain."""
        running = self.submit(PhBaseWorkChain, **self.ctx.inputs)
        self.report(f'running in serial, launching PhBaseWorkChain<{running.pk}>')
        return ToContext(workchain=running)

    def results(self):
        """Attach results to the workchain."""
        retrieved = self.ctx.workchain.out.retrieved
        self.out('retrieved', retrieved)
        self.report(f'workchain completed, output in {retrieved.__class__.__name__}<{retrieved.pk}>')
