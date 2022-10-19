# -*- coding: utf-8 -*-
"""Workchain to run Quantum ESPRESSO calculations that generate DoS and PDoS for a structure.

This requires four computations:

- SCF (pw.x), to generate the initial wavefunction.
- NSCF (pw.x), to generate eigenvalues, generally with a denser k-point mesh and tetrahedra occupations.
- Total DoS (dos.x), to generate total densities of state.
- Partial DoS (projwfc.x), to generate partial densities of state, by projecting wavefunctions onto atomic orbitals.

Additional functionality:

- Setting ``'align_to_fermi': True`` in the inputs will ensure that the energy range is centred around the Fermi
  energy when `Emin` and `Emax` are provided for both the `dos` and `projwfc` inputs. This is useful when you are only
  interested in a certain energy range around the Fermi energy. By default the energy range is extracted from the
  NSCF calculation.

Storage memory management:

The wavefunction file(s) created by the nscf calculation can get very large (>100Gb).
These files must be copied to dos and projwfc calculations, so storage memory limits can easily be exceeded.
If this is an issue, setting the input ``serial_clean`` to ``True`` will not run these calculations in parallel,
but instead run in serial and clean directories when they are no longer required:

- Run the scf workchain
- Run the nscf workchain, then clean the scf calculation directories
- Run the dos calculation, then clean its directory
- Run the projwfc calculation, then clean its directory

Setting the input ``clean_workdir`` to ``True``, will clean any remaining directories, after the whole workchain has
terminated.

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
from aiida import orm, plugins
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, if_
from aiida.orm.nodes.data.base import to_aiida_type
import jsonschema

from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.utils.mapping import prepare_process_inputs

from .protocols.utils import ProtocolMixin


def get_parameter_schema():
    """Return the ``PdosWorkChain`` input parameter schema."""
    return {
        '$schema': 'http://json-schema.org/draft-07/schema',
        'type': 'object',
        'required': ['DeltaE'],
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
                'minimum': 0
            },
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


def validate_inputs(value, _):
    """Validate the top level namespace.

    - Check that either the `scf` or `nscf.pw.parent_folder` inputs is provided.
    - Check that the `Emin`, `Emax` and `DeltaE` inputs are the same for the `dos` and `projwfc` namespaces.
    - Check that `Emin` and `Emax` are provided in case `align_to_fermi` is set to `True`.
    """
    # Check that either the `scf` input or `nscf.pw.parent_folder` is provided.
    import warnings
    if 'scf' in value and 'parent_folder' in value['nscf']['pw']:
        warnings.warn(
            'Both the `scf` and `nscf.pw.parent_folder` inputs were provided. The SCF calculation will '
            'be run with the inputs provided in `scf` and the `nscf.pw.parent_folder` will be ignored.'
        )
    elif not 'scf' in value and not 'parent_folder' in value['nscf']['pw']:
        return 'Specifying either the `scf` or `nscf.pw.parent_folder` input is required.'

    for par in ['Emin', 'Emax', 'DeltaE']:
        if value['dos']['parameters']['DOS'].get(par, None) != value['projwfc']['parameters']['PROJWFC'].get(par, None):
            return f'The `{par}`` parameter has to be equal for the `dos` and `projwfc` inputs.'

    if value.get('align_to_fermi', False):
        for par in ['Emin', 'Emax']:
            if value['dos']['parameters']['DOS'].get(par, None) is None:
                return f'The `{par}`` parameter must be set in case `align_to_fermi` is set to `True`.'


def validate_scf(value, _):
    """Validate the scf parameters."""
    parameters = value['pw']['parameters'].get_dict()
    if parameters.get('CONTROL', {}).get('calculation', 'scf') != 'scf':
        return '`CONTOL.calculation` in `scf.pw.parameters` is not set to `scf`.'


def validate_nscf(value, _):
    """Validate the nscf parameters."""
    parameters = value['pw']['parameters'].get_dict()
    if parameters.get('CONTROL', {}).get('calculation', 'scf') != 'nscf':
        return '`CONTOL.calculation` in `nscf.pw.parameters` is not set to `nscf`.'
    if parameters.get('SYSTEM', {}).get('occupations', None) != 'tetrahedra':
        return '`SYSTEM.occupations` in `nscf.pw.parameters` is not set to `tetrahedra`.'


def validate_dos(value, _):
    """Validate DOS parameters.

    - shared: Emin | Emax | DeltaE
    - dos.x only: ngauss | degauss | bz_sum
    - projwfc.x only: ngauss | degauss | pawproj | n_proj_boxes | irmin(3,n_proj_boxes) | irmax(3,n_proj_boxes)

    """
    jsonschema.validate(value['parameters'].get_dict()['DOS'], get_parameter_schema())


def validate_projwfc(value, _):
    """Validate DOS parameters.

    - shared: Emin | Emax | DeltaE
    - dos.x only: ngauss | degauss | bz_sum
    - projwfc.x only: ngauss | degauss | pawproj | n_proj_boxes | irmin(3,n_proj_boxes) | irmax(3,n_proj_boxes)

    """
    jsonschema.validate(value['parameters'].get_dict()['PROJWFC'], get_parameter_schema())


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


class PdosWorkChain(ProtocolMixin, WorkChain):
    """A WorkChain to compute Total & Partial Density of States of a structure, using Quantum Espresso."""

    @classmethod
    def define(cls, spec):
        # yapf: disable
        """Define the process specification."""
        super().define(spec)
        spec.input('structure', valid_type=orm.StructureData, help='The input structure.')
        spec.input(
            'serial_clean',
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
        spec.input(
            'align_to_fermi',
            valid_type=orm.Bool,
            serializer=to_aiida_type,
            default=lambda: orm.Bool(False),
            help=(
                'If true, Emin=>Emin-Efermi & Emax=>Emax-Efermi, where Efermi is taken from the `nscf` calculation. '
                'Note that it only makes sense to align `Emax` and `Emin` to the fermi level in case they are actually '
                'provided by in the `dos` and `projwfc` inputs, since otherwise the '
            )
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
            PwBaseWorkChain,
            namespace='nscf',
            exclude=('clean_workdir', 'pw.structure'),
            namespace_options={
                'help': 'Inputs for the `PwBaseWorkChain` of the `nscf` calculation.',
                'validator': validate_nscf
            }
        )
        spec.inputs['nscf']['pw'].validator = PwCalculation.validate_inputs_base
        spec.expose_inputs(
            DosCalculation,
            namespace='dos',
            exclude=('parent_folder',),
            namespace_options={
                'help': ('Input parameters for the `dos.x` calculation. Note that the `Emin`, `Emax` and `DeltaE` '
                         'values have to match with those in the `projwfc` inputs.'),
                'validator': validate_dos
            }
        )
        spec.expose_inputs(
            ProjwfcCalculation,
            namespace='projwfc',
            exclude=('parent_folder',),
            namespace_options={
                'help': ('Input parameters for the `projwfc.x` calculation. Note that the `Emin`, `Emax` and `DeltaE` '
                         'values have to match with those in the `dos` inputs.'),
                'validator': validate_projwfc
            }
        )
        spec.inputs.validator = validate_inputs

        spec.outline(
            cls.setup,
            if_(cls.should_run_scf)(
                cls.run_scf,
                cls.inspect_scf,
            ),
            cls.run_nscf,
            cls.inspect_nscf,
            if_(cls.serial_clean)(
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

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        from . import protocols
        return files(protocols) / 'pdos.yaml'

    @classmethod
    def get_builder_from_protocol(
        cls, pw_code, dos_code, projwfc_code, structure, protocol=None, overrides=None, options=None, **kwargs
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param pw_code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param dos_code: the ``Code`` instance configured for the ``quantumespresso.dos`` plugin.
        :param projwfc_code: the ``Code`` instance configured for the ``quantumespresso.projwfc`` plugin.
        :param structure: the ``StructureData`` instance to use.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param options: A dictionary of options that will be recursively set for the ``metadata.options`` input of all
            the ``CalcJobs`` that are nested in this work chain.
        :param kwargs: additional keyword arguments that will be passed to the ``get_builder_from_protocol`` of all the
            sub processes that are called by this workchain.
        :return: a process builder instance with all inputs defined ready for launch.
        """
        from aiida_quantumespresso.workflows.protocols.utils import recursive_merge

        inputs = cls.get_protocol_inputs(protocol, overrides)

        args = (pw_code, structure, protocol)
        scf = PwBaseWorkChain.get_builder_from_protocol(
            *args, overrides=inputs.get('scf', None), options=options, **kwargs
        )
        scf['pw'].pop('structure', None)
        scf.pop('clean_workdir', None)
        nscf = PwBaseWorkChain.get_builder_from_protocol(
            *args, overrides=inputs.get('nscf', None), options=options, **kwargs
        )
        nscf['pw'].pop('structure', None)
        nscf['pw']['parameters']['SYSTEM'].pop('smearing', None)
        nscf['pw']['parameters']['SYSTEM'].pop('degauss', None)
        nscf.pop('clean_workdir', None)

        metadata_dos = inputs.get('dos', {}).get('metadata', {'options': {}})
        metadata_projwfc = inputs.get('projwfc', {}).get('metadata', {'options': {}})

        if options:
            metadata_dos['options'] = recursive_merge(metadata_dos['options'], options)
            metadata_projwfc['options'] = recursive_merge(metadata_projwfc['options'], options)

        builder = cls.get_builder()
        builder.structure = structure
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])
        builder.scf = scf
        builder.nscf = nscf
        builder.dos.code = dos_code  # pylint: disable=no-member
        builder.dos.parameters = orm.Dict(inputs.get('dos', {}).get('parameters')) # pylint: disable=no-member
        builder.dos.metadata = metadata_dos  # pylint: disable=no-member
        builder.projwfc.code = projwfc_code  # pylint: disable=no-member
        builder.projwfc.parameters = orm.Dict(inputs.get('projwfc', {}).get('parameters')) # pylint: disable=no-member
        builder.projwfc.metadata = metadata_projwfc  # pylint: disable=no-member

        return builder

    def setup(self):
        """Initialize context variables that are used during the logical flow of the workchain."""
        self.ctx.serial_clean = 'serial_clean' in self.inputs and self.inputs.serial_clean.value
        self.ctx.dry_run = 'dry_run' in self.inputs and self.inputs.dry_run.value

    def serial_clean(self):
        """Return whether dos and projwfc calculations should be run in serial.

        The calculation remote folders will be cleaned before the next process step.
        """
        return self.ctx.serial_clean

    def should_run_scf(self):
        """Return whether the work chain should run an SCF calculation."""
        return 'scf' in self.inputs

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

    def run_nscf(self):
        """Run an NSCF calculation, to generate eigenvalues with a denser k-point mesh.

        This calculation modifies the base scf calculation inputs by:

        - Using the parent folder from the scf calculation.
        - Replacing the kpoints, if an alternative is specified for nscf.
        - Changing ``SYSTEM.occupations`` to 'tetrahedra'.
        - Changing ``SYSTEM.nosym`` to True, to avoid generation of additional k-points in low symmetry cases.
        - Replace the ``pw.metadata.options``, if an alternative is specified for nscf.

        """
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, 'nscf'))
        if 'scf' in self.inputs:
            inputs.pw.parent_folder = self.ctx.scf_parent_folder
        inputs.pw.structure = self.inputs.structure

        inputs.metadata.call_link_label = 'nscf'
        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)

        if self.ctx.dry_run:
            return inputs

        future = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launching NSCF PwBaseWorkChain<{future.pk}>')

        return ToContext(workchain_nscf=future)

    def inspect_nscf(self):
        """Verify that the NSCF calculation finished successfully."""
        workchain = self.ctx.workchain_nscf
        if not workchain.is_finished_ok:
            self.report(f'NSCF PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_NSCF

        if self.ctx.serial_clean:
            # we no longer require the scf remote folder, so can clean it
            cleaned_calcs = clean_workchain_calcs(self.ctx.workchain_scf)
            if cleaned_calcs:
                self.report(f"cleaned remote folders of SCF calculations: {' '.join(map(str, cleaned_calcs))}")

        self.ctx.nscf_emin = workchain.outputs.output_band.get_array('bands').min()
        self.ctx.nscf_emax = workchain.outputs.output_band.get_array('bands').max()
        self.ctx.nscf_parent_folder = workchain.outputs.remote_folder
        self.ctx.nscf_fermi = workchain.outputs.output_parameters.dict.fermi_energy

    def _generate_dos_inputs(self):
        """Run DOS calculation, to generate total Densities of State."""
        dos_inputs = AttributeDict(self.exposed_inputs(DosCalculation, 'dos'))
        dos_inputs.parent_folder = self.ctx.nscf_parent_folder
        dos_parameters = self.inputs.dos.parameters.get_dict()

        if dos_parameters.pop('align_to_fermi', False):
            dos_parameters['DOS']['Emin'] = dos_parameters['Emin'] + self.ctx.nscf_fermi
            dos_parameters['DOS']['Emax'] = dos_parameters['Emax'] + self.ctx.nscf_fermi

        dos_inputs.parameters = orm.Dict(dos_parameters)
        dos_inputs['metadata']['call_link_label'] = 'dos'
        return dos_inputs

    def _generate_projwfc_inputs(self):
        """Run Projwfc calculation, to generate partial Densities of State."""
        projwfc_inputs = AttributeDict(self.exposed_inputs(ProjwfcCalculation, 'projwfc'))
        projwfc_inputs.parent_folder = self.ctx.nscf_parent_folder
        projwfc_parameters = self.inputs.projwfc.parameters.get_dict()

        if projwfc_parameters.pop('align_to_fermi', False):
            projwfc_parameters['PROJWFC']['Emin'] = projwfc_parameters['Emin'] + self.ctx.nscf_fermi
            projwfc_parameters['PROJWFC']['Emax'] = projwfc_parameters['Emax'] + self.ctx.nscf_fermi

        projwfc_inputs.parameters = orm.Dict(projwfc_parameters)
        projwfc_inputs['metadata']['call_link_label'] = 'projwfc'
        return projwfc_inputs

    def run_dos_serial(self):
        """Run DOS calculation."""
        dos_inputs = self._generate_dos_inputs()

        if self.ctx.dry_run:
            return dos_inputs

        future_dos = self.submit(DosCalculation, **dos_inputs)
        self.report(f'launching DosCalculation<{future_dos.pk}>')
        return ToContext(calc_dos=future_dos)

    def inspect_dos_serial(self):
        """Verify that the DOS calculation finished successfully, then clean its remote directory."""
        calculation = self.ctx.calc_dos

        if not calculation.is_finished_ok:
            self.report(f'DosCalculation failed with exit status {calculation.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_DOS

        if self.ctx.serial_clean:
            # we no longer require the dos remote folder, so can clean it
            if clean_calcjob_remote(calculation):
                self.report(f'cleaned remote folder of DosCalculation<{calculation.pk}>')

    def run_projwfc_serial(self):
        """Run Projwfc calculation."""
        projwfc_inputs = self._generate_projwfc_inputs()

        if self.ctx.dry_run:
            return projwfc_inputs

        future_projwfc = self.submit(ProjwfcCalculation, **projwfc_inputs)
        self.report(f'launching ProjwfcCalculation<{future_projwfc.pk}>')
        return ToContext(calc_projwfc=future_projwfc)

    def inspect_projwfc_serial(self):
        """Verify that the Projwfc calculation finished successfully, then clean its remote directory."""
        calculation = self.ctx.calc_projwfc
        if not calculation.is_finished_ok:
            self.report(f'ProjwfcCalculation failed with exit status {calculation.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_PROJWFC

        if self.ctx.serial_clean:
            # we no longer require the projwfc remote folder, so can clean it
            if clean_calcjob_remote(calculation):
                self.report(f'cleaned remote folder of ProjwfcCalculation<{calculation.pk}>')

    def run_pdos_parallel(self):
        """Run DOS and Projwfc calculations in parallel."""
        dos_inputs = self._generate_dos_inputs()
        projwfc_inputs = self._generate_projwfc_inputs()

        if self.ctx.dry_run:
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

        self.out_many(self.exposed_outputs(self.ctx.workchain_nscf, PwBaseWorkChain, namespace='nscf'))
        self.out_many(self.exposed_outputs(self.ctx.calc_dos, DosCalculation, namespace='dos'))
        self.out_many(self.exposed_outputs(self.ctx.calc_projwfc, ProjwfcCalculation, namespace='projwfc'))

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super().on_terminated()

        if self.inputs.clean_workdir.value is False:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = clean_workchain_calcs(self.node)

        if cleaned_calcs:
            self.report(f"cleaned remote folders of calculations: {' '.join(map(str, cleaned_calcs))}")
