# -*- coding: utf-8 -*-
"""Workchain to automatically compute a band structure for a given structure using Quantum ESPRESSO pw.x."""
from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import WorkChain, ToContext
from aiida.plugins import WorkflowFactory

from aiida_quantumespresso.utils.protocols.pw import ProtocolManager
from aiida_quantumespresso.utils.pseudopotential import get_pseudos_from_dict
from aiida_quantumespresso.utils.resources import get_default_options

PwBandsWorkChain = WorkflowFactory('quantumespresso.pw.bands')


def validate_protocol(protocol_dict, ctx=None):  # pylint: disable=unused-argument
    """Check that the protocol is one for which we have a definition."""
    try:
        protocol_name = protocol_dict['name']
    except KeyError as exception:
        return 'Missing key `name` in protocol dictionary'
    try:
        ProtocolManager(protocol_name)
    except ValueError as exception:
        return str(exception)


class PwBandStructureWorkChain(WorkChain):
    """Workchain to automatically compute a band structure for a given structure using Quantum ESPRESSO pw.x."""

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('code', valid_type=orm.Code,
            help='The `pw.x` code to use for the `PwCalculations`.')
        spec.input('structure', valid_type=orm.StructureData,
            help='The input structure.')
        spec.input('options', valid_type=orm.Dict, required=False,
            help='Optional `options` to use for the `PwCalculations`.')
        spec.input('protocol', valid_type=orm.Dict, default=lambda: orm.Dict(dict={'name': 'theos-ht-1.0'}),
            help='The protocol to use for the workchain.', validator=validate_protocol)
        spec.expose_outputs(PwBandsWorkChain)
        spec.outline(
            cls.setup_protocol,
            cls.setup_parameters,
            cls.run_bands,
            cls.results,
        )
        spec.exit_code(201, 'ERROR_INVALID_INPUT_UNRECOGNIZED_KIND',
            message='Input `StructureData` contains an unsupported kind.')
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_BANDS',
            message='The `PwBandsWorkChain` sub process failed.')
        spec.output('primitive_structure', valid_type=orm.StructureData)
        spec.output('seekpath_parameters', valid_type=orm.Dict)
        spec.output('scf_parameters', valid_type=orm.Dict)
        spec.output('band_parameters', valid_type=orm.Dict)
        spec.output('band_structure', valid_type=orm.BandsData)
        # yapf: enable

    def _get_protocol(self):
        """Return a `ProtocolManager` instance and a dictionary of modifiers."""
        protocol_data = self.inputs.protocol.get_dict()
        protocol_name = protocol_data['name']
        protocol = ProtocolManager(protocol_name)

        protocol_modifiers = protocol_data.get('modifiers', {})

        return protocol, protocol_modifiers

    def setup_protocol(self):
        """Set up context variables and inputs for the `PwBandsWorkChain`.

        Based on the specified protocol, we define values for variables that affect the execution of the calculations.
        """
        protocol, protocol_modifiers = self._get_protocol()
        self.report(f'running the workchain with the "{protocol.name}" protocol')
        self.ctx.protocol = protocol.get_protocol_data(modifiers=protocol_modifiers)

    def setup_parameters(self):
        """Set up the default input parameters required for the `PwBandsWorkChain`."""
        ecutwfc = []
        ecutrho = []

        for kind in self.inputs.structure.get_kind_names():
            try:
                dual = self.ctx.protocol['pseudo_data'][kind]['dual']
                cutoff = self.ctx.protocol['pseudo_data'][kind]['cutoff']
                cutrho = dual * cutoff
                ecutwfc.append(cutoff)
                ecutrho.append(cutrho)
            except KeyError:
                self.report(f'failed to retrieve the cutoff or dual factor for {kind}')
                return self.exit_codes.ERROR_INVALID_INPUT_UNRECOGNIZED_KIND

        self.ctx.parameters = orm.Dict(
            dict={
                'CONTROL': {
                    'restart_mode': 'from_scratch',
                    'tstress': self.ctx.protocol['tstress'],
                    'tprnfor': self.ctx.protocol['tprnfor'],
                },
                'SYSTEM': {
                    'ecutwfc': max(ecutwfc),
                    'ecutrho': max(ecutrho),
                    'smearing': self.ctx.protocol['smearing'],
                    'degauss': self.ctx.protocol['degauss'],
                    'occupations': self.ctx.protocol['occupations'],
                },
                'ELECTRONS': {
                    'conv_thr': self.ctx.protocol['convergence_threshold_per_atom'] * len(self.inputs.structure.sites),
                }
            }
        )

    def run_bands(self):
        """Run the `PwBandsWorkChain` to compute the band structure."""

        def get_common_inputs():
            """Return the dictionary of inputs to be used as the basis for each `PwBaseWorkChain`."""
            protocol, protocol_modifiers = self._get_protocol()
            checked_pseudos = protocol.check_pseudos(
                modifier_name=protocol_modifiers.get('pseudo', None),
                pseudo_data=protocol_modifiers.get('pseudo_data', None)
            )
            known_pseudos = checked_pseudos['found']

            inputs = AttributeDict({
                'pw': {
                    'code': self.inputs.code,
                    'pseudos': get_pseudos_from_dict(self.inputs.structure, known_pseudos),
                    'parameters': self.ctx.parameters,
                    'metadata': {},
                }
            })

            if 'options' in self.inputs:
                inputs.pw.metadata.options = self.inputs.options.get_dict()
            else:
                inputs.pw.metadata.options = get_default_options(with_mpi=True)

            return inputs

        inputs = AttributeDict({
            'structure': self.inputs.structure,
            'relax': {
                'base': get_common_inputs(),
                'relaxation_scheme': orm.Str('vc-relax'),
                'meta_convergence': orm.Bool(self.ctx.protocol['meta_convergence']),
                'volume_convergence': orm.Float(self.ctx.protocol['volume_convergence']),
            },
            'scf': get_common_inputs(),
            'bands': get_common_inputs(),
        })

        inputs.relax.base.kpoints_distance = orm.Float(self.ctx.protocol['kpoints_mesh_density'])
        inputs.scf.kpoints_distance = orm.Float(self.ctx.protocol['kpoints_mesh_density'])
        inputs.bands.kpoints_distance = orm.Float(self.ctx.protocol['kpoints_distance_for_bands'])

        num_bands_factor = self.ctx.protocol.get('num_bands_factor', None)
        if num_bands_factor is not None:
            inputs.nbands_factor = orm.Float(num_bands_factor)

        running = self.submit(PwBandsWorkChain, **inputs)

        self.report(f'launching PwBandsWorkChain<{running.pk}>')

        return ToContext(workchain_bands=running)

    def results(self):
        """Attach the relevant output nodes from the band calculation to the workchain outputs for convenience."""
        workchain = self.ctx.workchain_bands

        if not self.ctx.workchain_bands.is_finished_ok:
            self.report(f'sub process PwBandsWorkChain<{workchain.pk}> failed')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_BANDS

        self.report('workchain successfully completed')
        link_labels = [
            'primitive_structure', 'seekpath_parameters', 'scf_parameters', 'band_parameters', 'band_structure'
        ]

        for link_triple in workchain.get_outgoing().all():
            if link_triple.link_label in link_labels:
                self.out(link_triple.link_label, link_triple.node)
