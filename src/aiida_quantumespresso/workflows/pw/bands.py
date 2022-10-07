# -*- coding: utf-8 -*-
"""Workchain to compute a band structure for a given structure using Quantum ESPRESSO pw.x."""
from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, if_
from aiida.plugins import WorkflowFactory

from aiida_quantumespresso.calculations.functions.seekpath_structure_analysis import seekpath_structure_analysis
from aiida_quantumespresso.utils.mapping import prepare_process_inputs

from ..protocols.utils import ProtocolMixin

PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
PwRelaxWorkChain = WorkflowFactory('quantumespresso.pw.relax')


def validate_inputs(inputs, ctx=None):  # pylint: disable=unused-argument
    """Validate the inputs of the entire input namespace."""
    # pylint: disable=no-member

    if 'nbands_factor' in inputs and 'nbnd' in inputs['bands']['pw']['parameters'].base.attributes.get('SYSTEM', {}):
        return PwBandsWorkChain.exit_codes.ERROR_INVALID_INPUT_NUMBER_OF_BANDS.message

    # Cannot specify both `bands_kpoints` and `bands_kpoints_distance`
    if all(key in inputs for key in ['bands_kpoints', 'bands_kpoints_distance']):
        return PwBandsWorkChain.exit_codes.ERROR_INVALID_INPUT_KPOINTS.message


class PwBandsWorkChain(ProtocolMixin, WorkChain):
    """Workchain to compute a band structure for a given structure using Quantum ESPRESSO pw.x.

    The logic for the computation of various parameters for the BANDS step is as follows:

    Number of bands:
        One can specify the number of bands to be used in the BANDS step either directly through the input parameters
        `bands.pw.parameters.SYSTEM.nbnd` or through `nbands_factor`. Note that specifying both is not allowed. When
        neither is specified nothing will be set by the work chain and the default of Quantum ESPRESSO will end up being
        used. If the `nbands_factor` is specified the maximum value of the following values will be used:

        * `nbnd` of the preceding SCF calculation
        * 0.5 * nspin * nelectrons * nbands_factor
        * 0.5 * nspin * nelectrons + 4 * nspin

    Kpoints:
        There are three options; specify either an existing `KpointsData` through `bands_kpoints`, or specify the
        `bands_kpoint_distance`, or specify neither. For the former those exact kpoints will be used for the BANDS step.
        In the two other cases, the structure will first be normalized using SeekPath and the path along high-symmetry
        k-points will be generated on that structure. The distance between kpoints for the path will be equal to that
        of `bands_kpoints_distance` or the SeekPath default if not specified.
    """

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.expose_inputs(PwRelaxWorkChain, namespace='relax', exclude=('clean_workdir', 'structure'),
            namespace_options={'required': False, 'populate_defaults': False,
            'help': 'Inputs for the `PwRelaxWorkChain`, if not specified at all, the relaxation step is skipped.'})
        spec.expose_inputs(PwBaseWorkChain, namespace='scf',
            exclude=('clean_workdir', 'pw.structure'),
            namespace_options={'help': 'Inputs for the `PwBaseWorkChain` for the SCF calculation.'})
        spec.expose_inputs(PwBaseWorkChain, namespace='bands',
            exclude=('clean_workdir', 'pw.structure', 'pw.kpoints', 'pw.kpoints_distance'),
            namespace_options={'help': 'Inputs for the `PwBaseWorkChain` for the BANDS calculation.'})
        spec.input('structure', valid_type=orm.StructureData, help='The inputs structure.')
        spec.input('clean_workdir', valid_type=orm.Bool, default=lambda: orm.Bool(False),
            help='If `True`, work directories of all called calculation will be cleaned at the end of execution.')
        spec.input('nbands_factor', valid_type=orm.Float, required=False,
            help='The number of bands for the BANDS calculation is that used for the SCF multiplied by this factor.')
        spec.input('bands_kpoints', valid_type=orm.KpointsData, required=False,
            help='Explicit kpoints to use for the BANDS calculation. Specify either this or `bands_kpoints_distance`.')
        spec.input('bands_kpoints_distance', valid_type=orm.Float, required=False,
            help='Minimum kpoints distance for the BANDS calculation. Specify either this or `bands_kpoints`.')
        spec.inputs.validator = validate_inputs
        spec.outline(
            cls.setup,
            if_(cls.should_run_relax)(
                cls.run_relax,
                cls.inspect_relax,
            ),
            if_(cls.should_run_seekpath)(
                cls.run_seekpath,
            ),
            cls.run_scf,
            cls.inspect_scf,
            cls.run_bands,
            cls.inspect_bands,
            cls.results,
        )
        spec.exit_code(201, 'ERROR_INVALID_INPUT_NUMBER_OF_BANDS',
            message='Cannot specify both `nbands_factor` and `bands.pw.parameters.SYSTEM.nbnd`.')
        spec.exit_code(202, 'ERROR_INVALID_INPUT_KPOINTS',
            message='Cannot specify both `bands_kpoints` and `bands_kpoints_distance`.')
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_RELAX',
            message='The PwRelaxWorkChain sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_SCF',
            message='The scf PwBasexWorkChain sub process failed')
        spec.exit_code(403, 'ERROR_SUB_PROCESS_FAILED_BANDS',
            message='The bands PwBasexWorkChain sub process failed')
        spec.output('primitive_structure', valid_type=orm.StructureData,
            required=False,
            help='The normalized and primitivized structure for which the bands are computed.')
        spec.output('seekpath_parameters', valid_type=orm.Dict,
            required=False,
            help='The parameters used in the SeeKpath call to normalize the input or relaxed structure.')
        spec.output('scf_parameters', valid_type=orm.Dict,
            help='The output parameters of the SCF `PwBaseWorkChain`.')
        spec.output('band_parameters', valid_type=orm.Dict,
            help='The output parameters of the BANDS `PwBaseWorkChain`.')
        spec.output('band_structure', valid_type=orm.BandsData,
            help='The computed band structure.')
        # yapf: enable

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        from ..protocols import pw as pw_protocols
        return files(pw_protocols) / 'bands.yaml'

    @classmethod
    def get_builder_from_protocol(cls, code, structure, protocol=None, overrides=None, options=None, **kwargs):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param structure: the ``StructureData`` instance to use.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the protocol.
        :param options: A dictionary of options that will be recursively set for the ``metadata.options`` input of all
            the ``CalcJobs`` that are nested in this work chain.
        :param kwargs: additional keyword arguments that will be passed to the ``get_builder_from_protocol`` of all the
            sub processes that are called by this workchain.
        :return: a process builder instance with all inputs defined ready for launch.
        """
        inputs = cls.get_protocol_inputs(protocol, overrides)

        args = (code, structure, protocol)
        relax = PwRelaxWorkChain.get_builder_from_protocol(
            *args, overrides=inputs.get('relax', None), options=options, **kwargs
        )
        scf = PwBaseWorkChain.get_builder_from_protocol(
            *args, overrides=inputs.get('scf', None), options=options, **kwargs
        )
        bands = PwBaseWorkChain.get_builder_from_protocol(
            *args, overrides=inputs.get('bands', None), options=options, **kwargs
        )

        relax.pop('structure', None)
        relax.pop('clean_workdir', None)
        relax.pop('base_final_scf', None)
        scf['pw'].pop('structure', None)
        scf.pop('clean_workdir', None)
        bands['pw'].pop('structure', None)
        bands.pop('clean_workdir', None)
        bands.pop('kpoints_distance', None)
        bands.pop('kpoints_force_parity', None)

        builder = cls.get_builder()
        builder.structure = structure
        builder.relax = relax
        builder.scf = scf
        builder.bands = bands
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])
        builder.nbands_factor = orm.Float(inputs['nbands_factor'])
        builder.bands_kpoints_distance = orm.Float(inputs['bands_kpoints_distance'])

        return builder

    def setup(self):
        """Define the current structure in the context to be the input structure."""
        self.ctx.current_structure = self.inputs.structure
        self.ctx.current_number_of_bands = None
        self.ctx.bands_kpoints = self.inputs.get('bands_kpoints', None)

    def should_run_relax(self):
        """If the 'relax' input namespace was specified, we relax the input structure."""
        return 'relax' in self.inputs

    def should_run_seekpath(self):
        """Seekpath should only be run if the `bands_kpoints` input is not specified."""
        return 'bands_kpoints' not in self.inputs

    def run_relax(self):
        """Run the PwRelaxWorkChain to run a relax PwCalculation."""
        inputs = AttributeDict(self.exposed_inputs(PwRelaxWorkChain, namespace='relax'))
        inputs.metadata.call_link_label = 'relax'
        inputs.structure = self.ctx.current_structure

        running = self.submit(PwRelaxWorkChain, **inputs)

        self.report(f'launching PwRelaxWorkChain<{running.pk}>')

        return ToContext(workchain_relax=running)

    def inspect_relax(self):
        """Verify that the PwRelaxWorkChain finished successfully."""
        workchain = self.ctx.workchain_relax

        if not workchain.is_finished_ok:
            self.report(f'PwRelaxWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        self.ctx.current_structure = workchain.outputs.output_structure
        self.ctx.current_number_of_bands = workchain.outputs.output_parameters.base.attributes.get('number_of_bands')

    def run_seekpath(self):
        """Run the structure through SeeKpath to get the normalized structure and path along high-symmetry k-points .

        This is only called if the `bands_kpoints` input was not specified.
        """
        inputs = {
            'reference_distance': self.inputs.get('bands_kpoints_distance', None),
            'metadata': {
                'call_link_label': 'seekpath'
            }
        }
        result = seekpath_structure_analysis(self.ctx.current_structure, **inputs)
        self.ctx.current_structure = result['primitive_structure']
        self.ctx.bands_kpoints = result['explicit_kpoints']

        self.out('primitive_structure', result['primitive_structure'])
        self.out('seekpath_parameters', result['parameters'])

    def run_scf(self):
        """Run the PwBaseWorkChain in scf mode on the primitive cell of (optionally relaxed) input structure."""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        inputs.metadata.call_link_label = 'scf'
        inputs.pw.structure = self.ctx.current_structure
        inputs.pw.parameters = inputs.pw.parameters.get_dict()
        inputs.pw.parameters.setdefault('CONTROL', {})['calculation'] = 'scf'

        # Make sure to carry the number of bands from the relax workchain if it was run and it wasn't explicitly defined
        # in the inputs. One of the base workchains in the relax workchain may have changed the number automatically in
        #  the sanity checks on band occupations.
        if self.ctx.current_number_of_bands:
            inputs.pw.parameters.setdefault('SYSTEM', {}).setdefault('nbnd', self.ctx.current_number_of_bands)

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launching PwBaseWorkChain<{running.pk}> in scf mode')

        return ToContext(workchain_scf=running)

    def inspect_scf(self):
        """Verify that the PwBaseWorkChain for the scf run finished successfully."""
        workchain = self.ctx.workchain_scf

        if not workchain.is_finished_ok:
            self.report(f'scf PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF

        self.ctx.current_folder = workchain.outputs.remote_folder
        self.ctx.current_number_of_bands = workchain.outputs.output_parameters.base.attributes.get('number_of_bands')

    def run_bands(self):
        """Run the PwBaseWorkChain in bands mode along the path of high-symmetry determined by seekpath."""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='bands'))
        inputs.metadata.call_link_label = 'bands'
        inputs.kpoints = self.ctx.bands_kpoints
        inputs.pw.structure = self.ctx.current_structure
        inputs.pw.parent_folder = self.ctx.current_folder
        inputs.pw.parameters = inputs.pw.parameters.get_dict()
        inputs.pw.parameters.setdefault('CONTROL', {})
        inputs.pw.parameters.setdefault('SYSTEM', {})
        inputs.pw.parameters.setdefault('ELECTRONS', {})

        # The following flags always have to be set in the parameters, regardless of what caller specified in the inputs
        inputs.pw.parameters['CONTROL']['calculation'] = 'bands'

        # Only set the following parameters if not directly explicitly defined in the inputs
        inputs.pw.parameters['ELECTRONS'].setdefault('diagonalization', 'cg')
        inputs.pw.parameters['ELECTRONS'].setdefault('diago_full_acc', True)

        # If `nbands_factor` is defined in the inputs we set the `nbnd` parameter
        if 'nbands_factor' in self.inputs:
            factor = self.inputs.nbands_factor.value
            parameters = self.ctx.workchain_scf.outputs.output_parameters.get_dict()
            if int(parameters['number_of_spin_components']) > 1:
                nspin_factor = 2
            else:
                nspin_factor = 1
            nbands = int(parameters['number_of_bands'])
            nelectron = int(parameters['number_of_electrons'])
            nbnd = max(
                int(0.5 * nelectron * nspin_factor * factor),
                int(0.5 * nelectron * nspin_factor) + 4 * nspin_factor, nbands
            )
            inputs.pw.parameters['SYSTEM']['nbnd'] = nbnd

        # Otherwise set the current number of bands, unless explicitly set in the inputs
        else:
            inputs.pw.parameters['SYSTEM'].setdefault('nbnd', self.ctx.current_number_of_bands)

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launching PwBaseWorkChain<{running.pk}> in bands mode')

        return ToContext(workchain_bands=running)

    def inspect_bands(self):
        """Verify that the PwBaseWorkChain for the bands run finished successfully."""
        workchain = self.ctx.workchain_bands

        if not workchain.is_finished_ok:
            self.report(f'bands PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_BANDS

    def results(self):
        """Attach the desired output nodes directly as outputs of the workchain."""
        self.report('workchain succesfully completed')
        self.out('scf_parameters', self.ctx.workchain_scf.outputs.output_parameters)
        self.out('band_parameters', self.ctx.workchain_bands.outputs.output_parameters)
        self.out('band_structure', self.ctx.workchain_bands.outputs.output_band)

    def on_terminated(self):
        """Clean the working directories of all child calculations if `clean_workdir=True` in the inputs."""
        super().on_terminated()

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
            self.report(f"cleaned remote folders of calculations: {' '.join(map(str, cleaned_calcs))}")
