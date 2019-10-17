# -*- coding: utf-8 -*-
"""Workchain to run Quantum ESPRESSO calculations that generate DoS and PDoS for a structure.

This requires four computations:

- SCF (pw.x), to generate the initial wavefunction
- NSCF (pw.x), to generate eigenvalues with a denser k-point mesh
- Total DoS (dos.x), to generate total densities of state
- Partial DoS (projwfc.x), to generate partial densities of state by projecting wavefunctions onto atomic orbitals

Related Resources:

- https://www.quantum-espresso.org/Doc/pw_user_guide/node10.html
- https://blog.levilentz.com/density-of-states-calculation/
- http://indico.ictp.it/event/7921/session/320/contribution/1261/material/0/0.pdf

"""
from __future__ import absolute_import

import jsonschema
from six.moves import map

from plumpy import ProcessSpec

from aiida import engine, orm, plugins
from aiida.common import AttributeDict, LinkType
from aiida.orm.nodes.data.base import to_aiida_type

from aiida_quantumespresso.utils.mapping import prepare_process_inputs


def validate_dos_parameters(node):
    """Validate DOS parameters

    - shared: ngauss | degauss | Emin | Emax | DeltaE
    - dos.x only: bz_sum
    - projwfc.x only: pawproj | n_proj_boxes | irmin(3,n_proj_boxes) | irmax(3,n_proj_boxes)

    """
    jsonschema.validate(
        node.get_dict(), {
            '$schema': 'http://json-schema.org/draft-07/schema',
            'type': 'object',
            'required': ['Emin', 'Emax', 'DeltaE'],
            'additionalProperties': False,
            'properties': {
                'Emin': {
                    'description': 'min energy (eV) for DOS plot',
                    'type': 'number'
                },
                'Emax': {
                    'description': 'max energy (eV) for DOS plot',
                    'type': 'number'
                },
                'DeltaE': {
                    'description': 'energy grid step (eV)',
                    'type': 'number',
                },
                'ngauss': {
                    'description': 'Type of gaussian broadening.',
                    'type': 'integer',
                    'enum': [0, 1, -1, -99]
                },
                'degauss': {
                    'description': 'gaussian broadening, Ry (not eV!)',
                    'type': 'number'
                },
                'dos_only': {
                    'type': 'object'
                },
                'projwfc_only': {
                    'type': 'object'
                }
            }
        }
    )


PwBaseWorkChain = plugins.WorkflowFactory('quantumespresso.pw.base')
DosCalculation = plugins.CalculationFactory('quantumespresso.dos')
ProjwfcCalculation = plugins.CalculationFactory('quantumespresso.projwfc')


class PdosWorkChain(engine.WorkChain):
    """A WorkChain to compute Total & Partial Density of States of a structure, using Quantum Espresso."""

    @classmethod
    def define(cls, spec):
        # type: (ProcessSpec) -> None
        # yapf: disable
        """Define the process specification."""
        super(PdosWorkChain, cls).define(spec)

        # Base SCF & NSCF inputs
        spec.expose_inputs(
            PwBaseWorkChain, namespace='base',
            exclude=('clean_workdir', 'pw.parent_folder'),
            namespace_options={
                'help': 'Inputs for the `PwBaseWorkChain` to run scf and nscf calculations.'})
        # TODO optionally allow PwCalculation output parent_folder? # pylint: disable=fixme

        # optional NSCF overrides
        spec.expose_inputs(
            PwBaseWorkChain, namespace='nscf',
            include=('kpoints', 'kpoints_distance', 'kpoints_force_parity', 'pw.metadata.options'),
            namespace_options={
                'help': 'Optional inputs for nscf calculation, that override those set in the `base` namespace'})
        spec.input('nscf.pw.metadata.options', valid_type=dict, required=False)

        # DOS & PDOS
        spec.input(
            'parameters',
            valid_type=orm.Dict,
            serializer=to_aiida_type,
            validator=validate_dos_parameters,
            help='Combined input parameters for dos.x and projwfc.x'
        )
        spec.expose_inputs(DosCalculation, namespace='dos', exclude=('parent_folder', 'parameters'))
        spec.expose_inputs(ProjwfcCalculation, namespace='projwfc', exclude=('parent_folder', 'parameters'))

        # additional
        spec.input(
            'clean_workdir',
            valid_type=orm.Bool,
            serializer=to_aiida_type,
            default=orm.Bool(False),
            help='If ``True``, work directories of all called calculation will be cleaned at the end of execution.'
        )
        spec.input(
            'test_run',
            valid_type=orm.Bool,
            required=False,
            serializer=to_aiida_type,
            help='Terminate workchain steps before submitting calculations (test purposes only).')

        spec.outline(
            cls.validate,
            cls.run_scf,
            cls.inspect_scf,
            cls.run_nscf,
            cls.inspect_nscf,
            cls.run_dos,
            cls.inspect_dos,
            cls.results,
        )

        spec.exit_code(202, 'ERROR_INVALID_INPUT_KPOINTS',
            message='Neither the `kpoints` nor the `kpoints_distance` input was specified for base or nscf namespaces.')
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_SCF',
            message='the SCF sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_NSCF',
            message='the NSCF sub process failed')
        spec.exit_code(403, 'ERROR_SUB_PROCESS_FAILED_DOS',
            message='the DOS sub process failed')
        spec.exit_code(404, 'ERROR_SUB_PROCESS_FAILED_PROJWFC',
            message='the PROJWFC sub process failed')

        spec.expose_outputs(DosCalculation, namespace='dos')
        spec.expose_outputs(ProjwfcCalculation, namespace='projwfc')

    def validate(self):
        """Initialize context variables that are used during the logical flow of the workchain."""

    def run_scf(self):
        """Run an SCF calculation, to generate the wavefunction."""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, 'base'))
        inputs.pw.parameters = inputs.pw.parameters.get_dict()
        inputs.pw.parameters.setdefault('CONTROL', {})
        inputs.pw.parameters['CONTROL']['calculation'] = 'scf'
        inputs.pw.parameters['CONTROL']['restart_mode'] = 'from_scratch'

        inputs.metadata.call_link_label = 'workchain_scf'
        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)

        if 'test_run' in self.inputs and self.inputs.test_run.value:
            return inputs

        future = self.submit(PwBaseWorkChain, **inputs)

        self.report('launching SCF PwBaseWorkChain<{}>'.format(future.pk))

        return engine.ToContext(workchain_scf=future)

    def inspect_scf(self):
        """Verify that the SCF calculation finished successfully."""
        workchain = self.ctx.workchain_scf
        if not workchain.is_finished_ok:
            self.report('SCF PwBaseWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF

        self.ctx.scf_parent_folder = workchain.outputs.remote_workdir

    def run_nscf(self):
        """Run an NSCF calculation, to generate eigenvalues with a denser k-point mesh.

        This calculation modifies the base scf calculation inputs by:

        - Using the parent folder from the scf calculation.
        - Replacing the kpoints, if an alternative is specified for nscf.
        - Changing ``SYSTEM.occupations`` to 'tetrahedra'.
        - Changing ``SYSTEM.nosym`` to True, to avoid generation of additional k-points in low symmetry cases.
        - Replace the ``pw.metadata.options``, if an alternative is specified for nscf.

        """
        nscf_inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, 'nscf'))
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, 'base'))
        inputs.pw.parent_folder = self.ctx.scf_parent_folder
        if 'kpoints' in nscf_inputs or 'kpoints_distance' in nscf_inputs:
            for key in ['kpoints', 'kpoints_distance', 'kpoints_force_parity']:
                inputs.pop(key, None)
                if key in nscf_inputs:
                    inputs[key] = nscf_inputs[key]
        inputs.pw.parameters = inputs.pw.parameters.get_dict()
        inputs.pw.parameters.setdefault('CONTROL', {})
        inputs.pw.parameters['CONTROL']['calculation'] = 'nscf'
        inputs.pw.parameters['CONTROL']['restart_mode'] = 'from_scratch'
        inputs.pw.parameters.setdefault('SYSTEM', {})
        inputs.pw.parameters['SYSTEM']['occupations'] = 'tetrahedra'
        inputs.pw.parameters['SYSTEM']['nosym'] = True
        try:
            inputs.pw.metadata.options = nscf_inputs.pw.metadata.options
        except AttributeError:
            pass

        inputs.metadata.call_link_label = 'workchain_scf'
        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)

        if 'test_run' in self.inputs and self.inputs.test_run.value:
            return inputs

        future = self.submit(PwBaseWorkChain, **inputs)

        self.report('launching NSCF PwBaseWorkChain<{}>'.format(future.pk))

        return engine.ToContext(workchain_nscf=future)

    def inspect_nscf(self):
        """Verify that the NSCF calculation finished successfully."""
        workchain = self.ctx.workchain_nscf
        if not workchain.is_finished_ok:
            self.report('NSCF PwBaseWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_NSCF

        self.ctx.nscf_parent_folder = workchain.outputs.remote_workdir

    def run_dos(self):
        """Run DOS and Projwfc calculations, to generate total/partial Densities of State."""
        dos_inputs = AttributeDict(self.exposed_inputs(DosCalculation, 'dos'))
        dos_inputs.parent_folder = self.ctx.nscf_parent_folder
        dos_dict = self.inputs.parameters.get_dict()
        dos_dict.pop('projwfc_only', None)
        for key, val in dos_dict.pop('dos_only', {}).items():
            if key not in dos_dict:
                dos_dict[key] = val
        dos_inputs.parameters = orm.Dict(dict={
            'DOS': dos_dict
        })
        dos_inputs['metadata']['call_link_label'] = 'calc_dos'

        projwfc_inputs = AttributeDict(self.exposed_inputs(ProjwfcCalculation, 'projwfc'))
        projwfc_inputs.parent_folder = self.ctx.nscf_parent_folder
        projwfc_dict = self.inputs.parameters.get_dict()
        projwfc_dict.pop('dos_only', None)
        for key, val in projwfc_dict.pop('projwfc_only', {}).items():
            if key not in projwfc_dict:
                projwfc_dict[key] = val
        projwfc_inputs.parameters = orm.Dict(dict={
            'PROJWFC': projwfc_dict
        })
        projwfc_inputs['metadata']['call_link_label'] = 'calc_projwfc'

        if 'test_run' in self.inputs and self.inputs.test_run.value:
            return dos_inputs, projwfc_inputs

        future_dos = self.submit(DosCalculation, **dos_inputs)

        self.report('launching DosCalculation<{}>'.format(future_dos.pk))

        self.to_context(**{'calc_dos': future_dos})

        future_projwfc = self.submit(ProjwfcCalculation, **projwfc_inputs)

        self.report('launching ProjwfcCalculation<{}>'.format(future_projwfc.pk))

        self.to_context(**{'calc_projwfc': future_projwfc})

    def inspect_dos(self):
        """Verify that the DOS and Projwfc calculations finished successfully."""
        calculation = self.ctx.calc_dos
        if not calculation.is_finished_ok:
            self.report('DosCalculation failed with exit status {}'.format(calculation.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_DOS

        calculation = self.ctx.calc_projwfc
        if not calculation.is_finished_ok:
            self.report('ProjwfcCalculation failed with exit status {}'.format(calculation.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_PROJWFC

    def results(self):
        """Attach the desired output nodes directly as outputs of the workchain."""
        self.report('workchain succesfully completed')
        # TODO exposed_outputs for CalcJobs is fixed in aiida-core v1.0.0b6 # pylint: disable=fixme
        # self.out_many(self.exposed_outputs(calc_node, process_class, namespace=pname))
        namespace_separator = self.spec().namespace_separator
        for link_triple in self.ctx.calc_dos.get_outgoing(link_type=LinkType.CREATE).link_triples:
            self.out('dos' + namespace_separator + link_triple.link_label, link_triple.node)
        for link_triple in self.ctx.calc_projwfc.get_outgoing(link_type=LinkType.CREATE).link_triples:
            self.out('projwfc' + namespace_separator + link_triple.link_label, link_triple.node)

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super(PdosWorkChain, self).on_terminated()

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
