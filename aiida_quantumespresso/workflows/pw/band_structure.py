# -*- coding: utf-8 -*-
import copy

from aiida.orm import Code
from aiida.orm.data.base import Str, Float, Bool
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm import WorkflowFactory
from aiida.work.workchain import WorkChain, ToContext

from aiida_quantumespresso.utils.mapping import update_mapping, prepare_process_inputs
from aiida_quantumespresso.utils.protocols.pw import ProtocolManager
from aiida_quantumespresso.utils.pseudopotential import get_pseudos_from_dict

PwBandsWorkChain = WorkflowFactory('quantumespresso.pw.bands')

          
class PwBandStructureWorkChain(WorkChain):
    """
    Workchain to relax and compute the band structure for a given input structure
    using Quantum ESPRESSO's pw.x
    """

    ERROR_INVALID_INPUT_UNRECOGNIZED_KIND = 1

    @classmethod
    def define(cls, spec):
        super(PwBandStructureWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        #spec.input('pseudo_family', valid_type=Str)
        spec.input('protocol', valid_type=ParameterData)
        spec.input('scf_options', valid_type=ParameterData)
        #spec.expose_inputs(PwBandsWorkChain, namespace='bands_workflow', include=('relax.base.options', 'scf.options', 'bands.options'))
        #del spec.inputs['bands_workflow']['bands']['code']
        #del spec.inputs['bands_workflow']['scf']['code']
        #del spec.inputs['bands_workflow']['relax']['base']['code']
        spec.outline(
            cls.setup_protocol,
            cls.setup_kpoints,
            cls.setup_parameters,
            cls.run_bands,
            cls.run_results,
        )
        
        spec.exit_code(101,'BANDS_WORKCHAIN_FAILED',
            message='The bands subworkchain failed')
 
        spec.output('primitive_structure', valid_type=StructureData)
        spec.output('seekpath_parameters', valid_type=ParameterData)
        spec.output('scf_parameters', valid_type=ParameterData)
        spec.output('band_parameters', valid_type=ParameterData)
        spec.output('band_structure', valid_type=BandsData)

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
        kpoints_mesh = KpointsData()
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
                return self.ERROR_INVALID_INPUT_UNRECOGNIZED_KIND

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

            return {
                'code': self.inputs.code,
                #'pseudo_family': self.inputs.pseudo_family,
                'pseudos': pseudos,
                'parameters': ParameterData(dict=self.ctx.parameters),
                'options': self.inputs.scf_options, # NOTE: for now we use these options for all three steps (relax, scf, bands)
                #'settings': {},
            } 

        relax_inputs = get_common_inputs()
        update_mapping(relax_inputs, {
            'kpoints': self.ctx.kpoints_mesh,
            'relaxation_scheme': Str('vc-relax'),  
            'meta_convergence': Bool(self.ctx.protocol['meta_convergence']),
            'volume_convergence': Float(self.ctx.protocol['volume_convergence']),
        })

        scf_inputs = get_common_inputs()
        update_mapping(scf_inputs, {
            'kpoints': self.ctx.kpoints_mesh,
        })

        #Updating the number of bands to compute
        num_bands_factor = self.ctx.protocol.get('num_bands_factor', None)
        bands_inputs = get_common_inputs()
        if num_bands_factor is not None:
            #Using a dict because it is converted to ParameterData in the next step
            bands_inputs['workchain_options'] = ParameterData(dict={'num_bands_factor':num_bands_factor})
        update_mapping(bands_inputs, {
            'kpoints_distance': Float(self.ctx.protocol['kpoints_distance_for_bands']),
        })

        # Final input preparation, wrapping dictionaries in ParameterData nodes
        inputs = {
            'structure': self.inputs.structure,
            'relax': {
                'base': relax_inputs,
            },
            'scf': scf_inputs,
            'bands': bands_inputs,
        }

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
                        'exit status: {})'.format(
                        wc_bands.process_state, wc_bands.exit_message, wc_bands.exit_status))
            exit_code = self.exit_codes.BANDS_WORKCHAIN_FAILED
            
        for link_label in ['primitive_structure', 'seekpath_parameters', 'scf_parameters', 'band_parameters', 'band_structure']:
            if link_label in self.ctx.workchain_bands.out:
                node = self.ctx.workchain_bands.get_outputs_dict()[link_label]
                self.out(link_label, node)
                self.report("attaching {}<{}> as an output node with label '{}'"
                    .format(node.__class__.__name__, node.pk, link_label))
        if exit_code:
            return exit_code
