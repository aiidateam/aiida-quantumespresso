# -*- coding: utf-8 -*-
"""Workchain to run Quantum ESPRESSO calculations that generate DoS and PDoS for a structure.

This requires four computations:

- SCF (pw.x), to generate the initial wavefunction.
- NSCF (pw.x), to generate eigenvalues, generally with a denser k-point mesh and tetrahedra occupations.
- Total DoS (dos.x), to generate total densities of state.
- Partial DoS (projwfc.x), to generate partial densities of state, by projecting wavefunctions onto atomic orbitals.

Additional functionality:

- Setting ``'align_to_fermi': True`` in the input ``parameters`` node,
  will ensure that the energy range is centred around the Fermi energy.

Storage memory management:

The wavefunction file(s) created by the nscf calculation can get very large (>100Gb).
These files must be copied to dos and projwfc calculations, so storage memory limits can easily be exceeded.
If this is an issue, setting the input ``serial_clean`` to ``True`` will not run these calculations in parallel,
but instead run in serial and clean directories when they are no longer required:

- Run the scf workchain
- Run the nscf workchain, then clean the scf calculation directories
- Run the dos calculation, then clean its directory
- Run the projwfc calculation, then clean its directory

Setting the input ``clean_workdir`` to ``True``, will clean any remaining directories,
after the whole workchain has terminated.

Also note that projwfc will fail if the scf/nscf calculations were run with a different number of procs/pools and
``wf_collect=.false.`` (this setting is deprecated in newer version of QE).

Related Resources:

- `Electronic structure calculations user guide <https://www.quantum-espresso.org/Doc/pw_user_guide/node10.html>`_
- `Density of States calculation blog <https://blog.levilentz.com/density-of-states-calculation/>`_
- `Quantum ESPRESSO tutorial slides <http://indico.ictp.it/event/7921/session/320/contribution/1261/material/0/0.pdf>`_

.. warning::

    For QE v6.1, there is an issue using ``tetrahedra`` occupations, as is recommended for ``nscf``,
    and both ``dos.x`` and ``projwfc.x`` will raise errors when reading the xml file
    (see `this post <https://lists.quantum-espresso.org/pipermail/users/2017-November/039656.html>`_).

"""
from __future__ import absolute_import

import jsonschema
from six.moves import map

from plumpy import ProcessSpec  # pylint: disable=unused-import

from aiida import engine, orm, plugins
from aiida.common import AttributeDict, LinkType
from aiida.orm.nodes.data.base import to_aiida_type

from aiida_quantumespresso.utils.mapping import prepare_process_inputs


def get_parameter_schema():
    """Return the ``PdosWorkChain`` input parameter schema."""
    return {
        '$schema': 'http://json-schema.org/draft-07/schema',
        'type': 'object',
        'definitions': {
            'only': {
                'type': 'object',
                'patternProperties': {
                    '^([Ee][mM][iI][nN]|[Ee][mM][aA][xX]|[Dd][Ee][lL][tT][aA][Ee])$': {
                        'description': 'Emin, Emax & DeltaE should not be specified here',
                        'not': {}
                    }
                },
                'properties': {
                    'ngauss': {
                        'description': 'Type of gaussian broadening.',
                        'type': 'integer',
                        'enum': [0, 1, -1, -99]
                    },
                    'degauss': {
                        'description': 'gaussian broadening, Ry (not eV!)',
                        'type': 'number',
                        'minimum': 0
                    },
                }
            }
        },
        'required': ['Emin', 'Emax', 'DeltaE'],
        'additionalProperties': False,
        'properties': {
            'align_to_fermi': {
                'description': 'if true, Emin=>Emin-Efermi & Emax=>Emax-Efermi (Efermi is taken from nscf)',
                'type': 'boolean'
            },
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
                'minimum': 0
            },
            'dos_only': {
                '$ref': '#/definitions/only'
            },
            'projwfc_only': {
                '$ref': '#/definitions/only'
            }
        }
    }


def validate_dos_parameters(node):
    """Validate DOS parameters.

    - shared: Emin | Emax | DeltaE
    - dos.x only: ngauss | degauss | bz_sum
    - projwfc.x only: ngauss | degauss | pawproj | n_proj_boxes | irmin(3,n_proj_boxes) | irmax(3,n_proj_boxes)

    """
    jsonschema.validate(node.get_dict(), get_parameter_schema())


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

        # Run options
        spec.input('serial_clean',
            valid_type=orm.Bool,
            serializer=to_aiida_type,
            required=False,
            help=('If ``True``, calculations will be run in serial, '
                  'and work directories will be cleaned before the next step.')
        )
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
            cls.setup,
            cls.run_scf,
            cls.inspect_scf,
            cls.run_nscf,
            cls.inspect_nscf,
            engine.if_(cls.clean_serial)(
                cls.run_dos_serial,
                cls.inspect_dos_serial,
                cls.run_projwfc_serial,
                cls.inspect_projwfc_serial
            ).else_(
                cls.run_pdos_parallel,
                cls.inspect_pdos_parallel,
            ),
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
        spec.exit_code(404, 'ERROR_SUB_PROCESS_FAILED_BOTH',
            message='both the DOS and PROJWFC sub process failed')

        spec.expose_outputs(PwBaseWorkChain, namespace='nscf')
        spec.expose_outputs(DosCalculation, namespace='dos')
        spec.expose_outputs(ProjwfcCalculation, namespace='projwfc')

    def setup(self):
        """Initialize context variables that are used during the logical flow of the workchain."""
        self.ctx.clean_serial = 'serial_clean' in self.inputs and self.inputs.serial_clean.value
        self.ctx.is_test_run = 'test_run' in self.inputs and self.inputs.test_run.value

    def clean_serial(self):
        """Return whether dos and projwfc calculations should be run in serial.

        The calculation remote folders will be cleaned before the next process step.
        """
        return self.ctx.clean_serial

    def run_scf(self):
        """Run an SCF calculation, to generate the wavefunction."""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, 'base'))
        inputs.pw.parameters = inputs.pw.parameters.get_dict()
        inputs.pw.parameters.setdefault('CONTROL', {})
        inputs.pw.parameters['CONTROL']['calculation'] = 'scf'
        inputs.pw.parameters['CONTROL']['restart_mode'] = 'from_scratch'

        inputs.metadata.call_link_label = 'workchain_scf'
        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)

        if self.ctx.is_test_run:
            return inputs

        future = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launching SCF PwBaseWorkChain<{future.pk}>')

        return engine.ToContext(workchain_scf=future)

    def inspect_scf(self):
        """Verify that the SCF calculation finished successfully."""
        workchain = self.ctx.workchain_scf
        if not workchain.is_finished_ok:
            self.report(f'SCF PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF

        self.ctx.scf_parent_folder = workchain.outputs.remote_folder

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
        # TODO have an input dict to control override of base parameters (with these as default)? # pylint: disable=fixme
        inputs.pw.parameters['SYSTEM']['occupations'] = 'tetrahedra'
        inputs.pw.parameters['SYSTEM']['nosym'] = True
        try:
            inputs.pw.metadata.options = nscf_inputs.pw.metadata.options
        except AttributeError:
            pass

        inputs.metadata.call_link_label = 'workchain_scf'
        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)

        if self.ctx.is_test_run:
            return inputs

        future = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launching NSCF PwBaseWorkChain<{future.pk}>')

        return engine.ToContext(workchain_nscf=future)

    def inspect_nscf(self):
        """Verify that the NSCF calculation finished successfully."""
        workchain = self.ctx.workchain_nscf
        if not workchain.is_finished_ok:
            self.report(f'NSCF PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_NSCF

        if self.ctx.clean_serial:
            # we no longer require the scf remote folder, so can clean it
            cleaned_calcs = clean_workchain_calcs(self.ctx.workchain_scf)
            if cleaned_calcs:
                self.report(f"cleaned remote folders of SCF calculations: {' '.join(map(str, cleaned_calcs))}")

        self.ctx.nscf_parent_folder = workchain.outputs.remote_folder
        self.ctx.nscf_fermi = workchain.outputs.output_parameters.dict.fermi_energy
        # TODO ensure fermi units are eV (and convert)?  # pylint: disable=fixme

    def _generate_dos_inputs(self):
        """Run DOS calculation, to generate total Densities of State."""
        dos_inputs = AttributeDict(self.exposed_inputs(DosCalculation, 'dos'))
        dos_inputs.parent_folder = self.ctx.nscf_parent_folder
        dos_dict = self.inputs.parameters.get_dict()
        dos_dict.pop('projwfc_only', None)
        if dos_dict.pop('align_to_fermi', False):
            dos_dict['Emin'] = dos_dict['Emin'] + self.ctx.nscf_fermi
            dos_dict['Emax'] = dos_dict['Emax'] + self.ctx.nscf_fermi
        for key, val in dos_dict.pop('dos_only', {}).items():
            dos_dict[key] = val
        dos_inputs.parameters = orm.Dict(dict={
            'DOS': dos_dict
        })
        dos_inputs['metadata']['call_link_label'] = 'calc_dos'
        return dos_inputs

    def _generate_projwfc_inputs(self):
        """Run Projwfc calculation, to generate partial Densities of State."""
        projwfc_inputs = AttributeDict(self.exposed_inputs(ProjwfcCalculation, 'projwfc'))
        projwfc_inputs.parent_folder = self.ctx.nscf_parent_folder
        projwfc_dict = self.inputs.parameters.get_dict()
        projwfc_dict.pop('dos_only', None)
        if projwfc_dict.pop('align_to_fermi', False):
            projwfc_dict['Emin'] = projwfc_dict['Emin'] + self.ctx.nscf_fermi
            projwfc_dict['Emax'] = projwfc_dict['Emax'] + self.ctx.nscf_fermi
        for key, val in projwfc_dict.pop('projwfc_only', {}).items():
            projwfc_dict[key] = val
        projwfc_inputs.parameters = orm.Dict(dict={
            'PROJWFC': projwfc_dict
        })
        projwfc_inputs['metadata']['call_link_label'] = 'calc_projwfc'
        return projwfc_inputs

    def run_dos_serial(self):
        """Run DOS calculation."""
        dos_inputs = self._generate_dos_inputs()
        if self.ctx.is_test_run:
            return dos_inputs
        future_dos = self.submit(DosCalculation, **dos_inputs)
        self.report(f'launching DosCalculation<{future_dos.pk}>')
        return engine.ToContext(calc_dos=future_dos)

    def inspect_dos_serial(self):
        """Verify that the DOS calculation finished successfully, then clean its remote directory."""
        calculation = self.ctx.calc_dos
        if not calculation.is_finished_ok:
            self.report(f'DosCalculation failed with exit status {calculation.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_DOS

        if self.ctx.clean_serial:
            # we no longer require the dos remote folder, so can clean it
            if clean_calcjob_remote(calculation):
                self.report(f'cleaned remote folder of DosCalculation<{calculation.pk}>')

    def run_projwfc_serial(self):
        """Run Projwfc calculation."""
        projwfc_inputs = self._generate_projwfc_inputs()
        if self.ctx.is_test_run:
            return projwfc_inputs
        future_projwfc = self.submit(ProjwfcCalculation, **projwfc_inputs)
        self.report(f'launching ProjwfcCalculation<{future_projwfc.pk}>')
        return engine.ToContext(calc_projwfc=future_projwfc)

    def inspect_projwfc_serial(self):
        """Verify that the Projwfc calculation finished successfully, then clean its remote directory."""
        calculation = self.ctx.calc_projwfc
        if not calculation.is_finished_ok:
            self.report(f'ProjwfcCalculation failed with exit status {calculation.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_PROJWFC

        if self.ctx.clean_serial:
            # we no longer require the projwfc remote folder, so can clean it
            if clean_calcjob_remote(calculation):
                self.report(f'cleaned remote folder of ProjwfcCalculation<{calculation.pk}>')

    def run_pdos_parallel(self):
        """Run DOS and Projwfc calculations in parallel."""
        dos_inputs = self._generate_dos_inputs()
        projwfc_inputs = self._generate_projwfc_inputs()

        if self.ctx.is_test_run:
            return dos_inputs, projwfc_inputs

        future_dos = self.submit(DosCalculation, **dos_inputs)
        self.report(f'launching DosCalculation<{future_dos.pk}>')
        self.to_context(**{'calc_dos': future_dos})

        future_projwfc = self.submit(ProjwfcCalculation, **projwfc_inputs)
        self.report(f'launching ProjwfcCalculation<{future_projwfc.pk}>')
        self.to_context(**{'calc_projwfc': future_projwfc})

    def inspect_pdos_parallel(self):
        """Verify that the DOS and Projwfc calculations finished successfully."""
        error_codes = []

        calculation = self.ctx.calc_dos
        if not calculation.is_finished_ok:
            self.report(f'DosCalculation failed with exit status {calculation.exit_status}')
            error_codes.append(self.exit_codes.ERROR_SUB_PROCESS_FAILED_DOS)

        calculation = self.ctx.calc_projwfc
        if not calculation.is_finished_ok:
            self.report(f'ProjwfcCalculation failed with exit status {calculation.exit_status}')
            error_codes.append(self.exit_codes.ERROR_SUB_PROCESS_FAILED_PROJWFC)

        if len(error_codes) > 1:
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_BOTH
        if len(error_codes) == 1:
            return error_codes[0]

    def results(self):
        """Attach the desired output nodes directly as outputs of the workchain."""
        self.report('workchain successfully completed')
        # TODO exposed_outputs for CalcJobs is fixed in aiida-core v1.0.0b6 # pylint: disable=fixme
        # self.out_many(self.exposed_outputs(calc_node, process_class, namespace=pname))
        namespace_separator = self.spec().namespace_separator
        for link_triple in self.ctx.workchain_nscf.get_outgoing(link_type=LinkType.RETURN).link_triples:
            self.out('nscf' + namespace_separator + link_triple.link_label, link_triple.node)
        for link_triple in self.ctx.calc_dos.get_outgoing(link_type=LinkType.CREATE).link_triples:
            self.out('dos' + namespace_separator + link_triple.link_label, link_triple.node)
        for link_triple in self.ctx.calc_projwfc.get_outgoing(link_type=LinkType.CREATE).link_triples:
            self.out('projwfc' + namespace_separator + link_triple.link_label, link_triple.node)

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super().on_terminated()

        if self.inputs.clean_workdir.value is False:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = clean_workchain_calcs(self.node)

        if cleaned_calcs:
            self.report(f"cleaned remote folders of calculations: {' '.join(map(str, cleaned_calcs))}")
