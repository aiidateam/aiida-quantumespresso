# -*- coding: utf-8 -*-
from __future__ import absolute_import

from aiida import orm
from aiida.engine import WorkChain, ToContext
from aiida.plugins import WorkflowFactory

from aiida_quantumespresso.utils.mapping import update_mapping
from aiida_quantumespresso.utils.protocols.pw import ProtocolManager
from aiida_quantumespresso.utils.pseudopotential import get_pseudos_from_dict
from aiida_quantumespresso.utils.resources import get_default_options

PwBandsWorkChain = WorkflowFactory('quantumespresso.pw.bands')


def _validate_protocol(protocol_dict):
    """ Check that the protocol is one for which we have a definition. """
    try:
        protocol_name = protocol_dict['name']
    except KeyError as e:
        return "Couldn't find key " + str(e) + " in protocol dictionary"
    try:
        protocol = ProtocolManager(protocol_name)
    except ValueError as e:
        return str(e)  # "Unknown protocol '{}'".format(name)


class PwBandStructureWorkChain(WorkChain):
    """Workchain to automatically compute a band structure for a given structure using Quantum ESPRESSO pw.x"""

    @classmethod
    def define(cls, spec):
        super(PwBandStructureWorkChain, cls).define(spec)
        spec.input('code', valid_type=orm.Code)
        spec.input('structure', valid_type=orm.StructureData)
        spec.input('protocol', valid_type=orm.Dict, default=orm.Dict(dict={'name':'theos-ht-1.0'}),
                   validator=_validate_protocol)
        spec.input('scf_options', valid_type=orm.Dict, default=orm.Dict(dict=get_default_options(with_mpi=True)))
        spec.outline(
            cls.setup_protocol,
            cls.setup_kpoints,
            cls.setup_parameters,
            cls.run_bands,
            cls.run_results,
        )
        spec.exit_code(101, 'ERROR_INVALID_INPUT_UNRECOGNIZED_KIND', message='The bands subworkchain failed')
        spec.exit_code(102, 'ERROR_SUB_PROCESS_FAILED_BANDS', message='the bands PwBasexWorkChain sub process failed')
        spec.output('primitive_structure', valid_type=orm.StructureData)
        spec.output('seekpath_parameters', valid_type=orm.Dict)
        spec.output('scf_parameters', valid_type=orm.Dict)
        spec.output('band_parameters', valid_type=orm.Dict)
        spec.output('band_structure', valid_type=orm.BandsData)

    def _get_protocol(self):
        """
        Return a ProtocolManager class and a dictionary of modifiers
        """
        protocol_data = self.inputs.protocol.get_dict()
        protocol_name = protocol_data['name']
        protocol = ProtocolManager(protocol_name)

        protocol_modifiers = protocol_data.get('modifiers', {})

        return protocol, protocol_modifiers

    def setup_protocol(self):
        """
        Setup of context variables and inputs for the PwBandsWorkChain. Based on the specified
        protocol, we define values for variables that affect the execution of the calculations
        """
        protocol, protocol_modifiers = self._get_protocol()
        self.report('running the workchain in the "{}" protocol'.format(protocol.name))
        self.ctx.protocol = protocol.get_protocol_data(modifiers=protocol_modifiers)

    def setup_kpoints(self):
        """
        Define the k-point mesh for the relax and scf calculations.
        """
        kpoints_mesh = orm.KpointsData()
        kpoints_mesh.set_cell_from_structure(self.inputs.structure)
        kpoints_mesh.set_kpoints_mesh_from_density(
            distance=self.ctx.protocol['kpoints_mesh_density'],
            offset=self.ctx.protocol['kpoints_mesh_offset']
        )

        self.ctx.kpoints_mesh = kpoints_mesh

    def setup_parameters(self):
        """
        Setup the default input parameters required for the PwBandsWorkChain
        """
        structure = self.inputs.structure
        ecutwfc = []
        ecutrho = []

        for kind in structure.get_kind_names():
            try:
                dual = self.ctx.protocol['pseudo_data'][kind]['dual']
                cutoff = self.ctx.protocol['pseudo_data'][kind]['cutoff']
                cutrho = dual * cutoff
                ecutwfc.append(cutoff)
                ecutrho.append(cutrho)
            except KeyError:
                self.report('failed to retrieve the cutoff or dual factor for {}'.format(kind))
                return self.exit_codes.ERROR_INVALID_INPUT_UNRECOGNIZED_KIND

        natoms = len(structure.sites)
        conv_thr = self.ctx.protocol['convergence_threshold_per_atom'] * natoms

        self.ctx.parameters = {
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
                'conv_thr': conv_thr,
            }
        }

    def run_bands(self):
        """
        Run the PwBandsWorkChain to compute the band structure
        """
        def get_common_inputs():

            protocol, protocol_modifiers = self._get_protocol()
            checked_pseudos = protocol.check_pseudos(
                modifier_name=protocol_modifiers.get('pseudo', None),
                pseudo_data=protocol_modifiers.get('pseudo_data', None))
            known_pseudos = checked_pseudos['found']

            pseudos = get_pseudos_from_dict(self.inputs.structure, known_pseudos)

            res = {
                'code': self.inputs.code,
                'pseudos': pseudos,
                'parameters': orm.Dict(dict=self.ctx.parameters),
            }

            if 'scf_options' in self.inputs:
                res['options'] = self.inputs.scf_options

            return res

        relax_inputs = {
            'base': get_common_inputs(),
            'relaxation_scheme': orm.Str('vc-relax'),
            'meta_convergence': orm.Bool(self.ctx.protocol['meta_convergence']),
            'volume_convergence': orm.Float(self.ctx.protocol['volume_convergence']),
        }
        relax_inputs['base']['kpoints'] = self.ctx.kpoints_mesh

        scf_inputs = get_common_inputs()
        update_mapping(scf_inputs, {
            'kpoints': self.ctx.kpoints_mesh,
        })

        bands_inputs = get_common_inputs()

        update_mapping(bands_inputs, {
            'kpoints_distance': orm.Float(self.ctx.protocol['kpoints_distance_for_bands']),
        })

        # Final input preparation, wrapping dictionaries in Dict nodes
        inputs = {
            'structure': self.inputs.structure,
            'relax': relax_inputs,
            'scf': scf_inputs,
            'bands': bands_inputs,
        }

        # Updating the number of bands to compute
        num_bands_factor = self.ctx.protocol.get('num_bands_factor', None)
        if num_bands_factor is not None:
            inputs['nbands_factor'] = orm.Float(num_bands_factor)

        running = self.submit(PwBandsWorkChain, **inputs)

        self.report('launching PwBandsWorkChain<{}>'.format(running.pk))

        return ToContext(workchain_bands=running)

    def run_results(self):
        """
        Attach the relevant output nodes from the band calculation to the workchain outputs
        for convenience
        """
        exit_code = 0

        if self.ctx.workchain_bands.is_finished_ok:
            self.report('workchain succesfully completed')
        else:
            wc_bands = self.ctx.workchain_bands
            self.report('bands sub-workchain failed (process state: {}, '
                        'exit message: {}, '
                        'exit status: {})'.format(wc_bands.process_state, wc_bands.exit_message, wc_bands.exit_status))
            exit_code = self.exit_codes.ERROR_SUB_PROCESS_FAILED_BANDS

        link_labels = [
            'primitive_structure',
            'seekpath_parameters',
            'scf_parameters',
            'band_parameters',
            'band_structure'
        ]

        for link_triple in self.ctx.workchain_bands.get_outgoing().all():
            if link_triple.link_label in link_labels:
                self.out(link_triple.link_label, link_triple.node)
                self.report("attaching {}<{}> as an output node with label '{}'"
                    .format(link_triple.node.__class__.__name__, link_triple.node.pk, link_triple.link_label))

        if exit_code:
            return exit_code
