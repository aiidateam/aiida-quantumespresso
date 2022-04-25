# -*- coding: utf-8 -*-
"""Workchain to run Quantum ESPRESSO calculations that generate phonon dispersion for a structure.

This requires four computations:

- SCF (pw.x), to generate the initial wavefunction.
- PH (pw.x), to generate dynamical matrices for a uniform mesh of q vectors.
- Q2R (q2r.x), to calculate of interatomic force constants.
- MAYDYN (matdyn.x), to generate Fourier interpolation for various q points.
"""
from aiida import orm, plugins
from aiida.common import AttributeDict
from aiida.engine import WorkChain, ToContext
from aiida.orm.nodes.data.base import to_aiida_type

from aiida_quantumespresso.utils.mapping import prepare_process_inputs

from .protocols.utils import ProtocolMixin

def validate_scf(value, _):
    """Validate the scf parameters."""
    parameters = value['pw']['parameters'].get_dict()
    if parameters.get('CONTROL', {}).get('calculation', 'scf') != 'scf':
        return '`CONTOL.calculation` in `scf.pw.parameters` is not set to `scf`.'

def clean_calcjob_remote(node):
    """Clean the remote directory of a ``CalcJobNode``."""
    cleaned = False
    try:
        node.outputs.remote_folder._clean()  # pylint: disable=protected-access
        cleaned = True
    except (IOError, OSError, KeyError):
        pass
    return cleaned

def clean_workchain_calcs(workchain):
    """Clean all remote directories of a workchain's descendant calculations."""
    cleaned_calcs = []

    for called_descendant in workchain.called_descendants:
        if isinstance(called_descendant, orm.CalcJobNode):
            if clean_calcjob_remote(called_descendant):
                cleaned_calcs.append(called_descendant.pk)

    return cleaned_calcs

PwBaseWorkChain = plugins.WorkflowFactory('quantumespresso.pw.base')
PhBaseWorkChain = plugins.WorkflowFactory('quantumespresso.ph.base')
Q2rBaseWorkChain = plugins.WorkflowFactory('quantumespresso.q2r.base')
MatdynBaseWorkChain = plugins.WorkflowFactory('quantumespresso.matdyn.base')

class PhononDispersionWorkChain(ProtocolMixin, WorkChain):
    """A WorkChain to compute phonon dispersion using Quamtum Espresso."""

    @classmethod
    def define(cls, spec):
        # yapf: disable
        """Define the process specification."""
        super().define(spec)
        spec.input('structure', valid_type=orm.StructureData, help='The input structure.')
        # spec.input(
        #     'serial_clean'
        # )
        spec.input(
            'clean_workdir',
            valid_type=orm.Bool,
            serializer=to_aiida_type,
            default=lambda: orm.Bool(False),
            help='If ``True``, work directories of all called calculation will be cleaned at the end of execution.'
        )
        spec.input(
            'dry_run',
            valid_type=orm.Bool,
            serializer=to_aiida_type,
            required=False,
            help='Terminate workchain steps before submitting calculations (test purposes only).'
        )
        spec.expose_inputs(
            PwBaseWorkChain,
            namespace='scf',
            exclude=('clean_workdir', 'pw.structure', 'pw.parent_folder'),
            namespace_options={
                'help': 'Inputs for the `PwBaseWorkChain` of the `scf` calculation.',
                'validator': validate_scf,
                'required': False,
                'populate_defaults': False,
            }
        )
        spec.expose_inputs(
            PhBaseWorkChain,
            namespace='ph',
            exclude=('clean_workdir', 'ph.parent_folder'),
            namespace_options={
                'help': 'Inputs for the `PhBaseWorkChain` of the `ph` calculation.',
            }
        )
        spec.expose_inputs(
            Q2rBaseWorkChain,
            namespace='q2r',
            exclude=('clean_workdir', 'q2r.parent_folder'),
            namespace_options={
                'help': 'Inputs for the `Q2rBaseWorkChain` of the `q2r` calculation.',
            }
        )
        spec.expose_inputs(
            MatdynBaseWorkChain,
            namespace='matdyn',
            exclude=('clean_workdir', 'matdyn.parent_folder'),
            namespace_options={
                'help': 'Inputs for the `MatdynBaseWorkChain` of the `matdyn` calculation.',
            }
        )

        spec.outline(
            cls.setup,
            cls.run_scf,
            cls.inspect_scf,
            cls.run_ph,
            cls.inspect_ph,
            cls.run_q2r,
            cls.inspect_q2r,
            cls.run_matdyn,
            cls.inspect_matdyn,
            cls.results,
        )

        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_SCF',
            message='the SCF sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_PH',
            message='the PH sub process failed')
        spec.exit_code(403, 'ERROR_SUB_PROCESS_FAILED_Q2R',
            message='the Q2R sub process failed')
        spec.exit_code(404, 'ERROR_SUB_PROCESS_FAILED_MATDYN',
            message='the MATDYN sub process failed')

        spec.expose_outputs(PhBaseWorkChain, namespace='ph')
        spec.expose_outputs(MatdynBaseWorkChain, namespace='matdyn')
        # yapf: enable

    def setup(self):
        """Initialize context variables that are used during the logical flow of the workchain."""
        self.ctx.dry_run = 'dry_run' in self.inputs and self.inputs.dry_run.value

    def run_scf(self):
        """Run an SCF calculation, to generate the wavefunction."""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, 'scf'))
        inputs.pw.structure = self.inputs.structure

        inputs.metadata.call_link_label = 'scf'
        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)

        if self.ctx.dry_run:
            return inputs

        future = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launching SCF PwBaseWorkChain<{future.pk}>')

        return ToContext(workchain_scf=future)

    def inspect_scf(self):
        """Verify that the SCF calculation finished successfully."""
        workchain = self.ctx.workchain_scf
        if not workchain.is_finished_ok:
            self.report(f'SCF PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF

        self.ctx.scf_parent_folder = workchain.outputs.remote_folder

    def run_ph(self):
        """Run a ph calculation, to generate dynamical matrices for a uniform mesh of q vectors."""
        inputs = AttributeDict(self.exposed_inputs(PhBaseWorkChain, 'ph'))
        inputs.ph.parent_folder = self.ctx.scf_paret_folder

        inputs.metadata.call_link_label = 'ph'
        inputs = prepare_process_inputs(PhBaseWorkChain, inputs)

        if self.ctx.dry_run:
            return inputs

        future = self.submit(PhBaseWorkChain, **inputs)

        self.report(f'launching PH PhBaseWorkChain<{future.pk}>')

        return ToContext(workchain_ph=future)

    def inspect_ph(self):
        """Verify that the PH calculation finished successfully."""
        workchain = self.ctx.workchain_ph
        if not workchain.is_finished_ok:
            self.report(f'PH PhBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_PH

        self.ctx.ph_parent_folder = workchain.outputs.remote_folder

    def run_q2r(self):
        pass

    def inspect_q2r(self):
        pass

    def run_matdyn(self):
        pass

    def inspect_matdyn(self):
        pass

    def results(self):
        """Attach the desired output nodes directly as outputs of the workchain."""
        self.report('workchain successfully completed')

        self.out_many(self.exposed_outputs(self.ctx.workchain_ph, PhBaseWorkChain, namespace='ph'))
        self.out_many(self.exposed_outputs(self.ctx.workchain_matdyn, MatdynBaseWorkChain, namespace='matdyn'))

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super().on_terminated()

        if self.inputs.clean_workdir.value is False:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = clean_workchain_calcs(self.node)

        if cleaned_calcs:
            self.report(f"cleaned remote folders of calculations: {' '.join(map(str, cleaned_calcs))}")
