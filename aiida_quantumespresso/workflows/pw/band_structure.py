# -*- coding: utf-8 -*-
from aiida.orm import Code
from aiida.orm.data.base import Str, Float, Bool
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.data.array.kpoints import KpointsData
from aiida.orm.utils import WorkflowFactory
from aiida.work.run import submit
from aiida.work.workchain import WorkChain, ToContext


PwBandsWorkChain = WorkflowFactory('quantumespresso.pw.bands')


class PwBandStructureWorkChain(WorkChain):
    """
    Workchain to relax and compute the band structure for a given input structure
    using Quantum ESPRESSO's pw.x
    """

    @classmethod
    def define(cls, spec):
        super(PwBandStructureWorkChain, cls).define(spec)
        spec.input('code', valid_type=Code)
        spec.input('structure', valid_type=StructureData)
        spec.input('pseudo_family', valid_type=Str)
        spec.input('protocol', valid_type=Str, default=Str('standard'))
        spec.outline(
            cls.setup_protocol,
            cls.setup_kpoints,
            cls.setup_parameters,
            cls.run_bands,
            cls.run_results,
        )
        spec.output('primitive_structure', valid_type=StructureData)
        spec.output('seekpath_parameters', valid_type=ParameterData)
        spec.output('scf_parameters', valid_type=ParameterData)
        spec.output('band_parameters', valid_type=ParameterData)
        spec.output('band_structure', valid_type=BandsData)

    def setup_protocol(self):
        """
        Setup of context variables and inputs for the PwBandsWorkChain. Based on the specified
        protocol, we define values for variables that affect the execution of the calculations
        """
        self.ctx.inputs = {
            'code': self.inputs.code,
            'structure': self.inputs.structure,
            'pseudo_family': self.inputs.pseudo_family,
            'parameters': {},
            'settings': {},
            'options': {
                'resources': {
                    'num_machines': 1
                },
                'max_wallclock_seconds': 1800,
            },
        }

        if self.inputs.protocol == 'standard':
            self.report('running the workchain in the "{}" protocol'.format(self.inputs.protocol.value))
            self.ctx.protocol = {
                'kpoints_mesh_offset': [0., 0., 0.],
                'kpoints_mesh_density': 0.2,
                'convergence_threshold': 2.E-06,
                'smearing': 'marzari-vanderbilt',
                'degauss': 0.02,
                'occupations': 'smearing',
                'tstress': True,
                'pseudo_data': {
                    'H':  {'cutoff': 55,  'dual': 8,  'pseudo': '031US'},
                    'He': {'cutoff': 55,  'dual': 4,  'pseudo': 'SG15'},
                    'Li': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.4'},
                    'Be': {'cutoff': 40,  'dual': 8,  'pseudo': 'GBRV-1.4'},
                    'B':  {'cutoff': 40,  'dual': 8,  'pseudo': '031PAW'},
                    'C':  {'cutoff': 50,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'N':  {'cutoff': 55,  'dual': 8,  'pseudo': 'THEOS'},
                    'O':  {'cutoff': 45,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'F':  {'cutoff': 50,  'dual': 8,  'pseudo': 'GBRV-1.4'},
                    'Ne': {'cutoff': 200, 'dual': 8,  'pseudo': '100PAW'},
                    'Na': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Mg': {'cutoff': 35,  'dual': 8,  'pseudo': 'GBRV-1.4'},
                    'Al': {'cutoff': 30,  'dual': 8,  'pseudo': '100PAW'},
                    'Si': {'cutoff': 30,  'dual': 8,  'pseudo': '100US'},
                    'P':  {'cutoff': 30,  'dual': 8,  'pseudo': '100US'},
                    'S':  {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Cl': {'cutoff': 35,  'dual': 8,  'pseudo': 'GBRV-1.4'},
                    'Ar': {'cutoff': 120, 'dual': 8,  'pseudo': '100US'},
                    'K':  {'cutoff': 50,  'dual': 8,  'pseudo': '100US'},
                    'Ca': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Sc': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Ti': {'cutoff': 35,  'dual': 8,  'pseudo': 'GBRV-1.4'},
                    'V':  {'cutoff': 40,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Cr': {'cutoff': 40,  'dual': 8,  'pseudo': 'GBRV-1.5'},
                    'Mn': {'cutoff': 70,  'dual': 12, 'pseudo': '031PAW'},
                    'Fe': {'cutoff': 90,  'dual': 12, 'pseudo': '031PAW'},
                    'Co': {'cutoff': 55,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Ni': {'cutoff': 45,  'dual': 8,  'pseudo': 'GBRV-1.4'},
                    'Cu': {'cutoff': 40,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Zn': {'cutoff': 40,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Ga': {'cutoff': 35,  'dual': 8,  'pseudo': '031US'},
                    'Ge': {'cutoff': 40,  'dual': 8,  'pseudo': '100PAW'},
                    'As': {'cutoff': 30,  'dual': 8,  'pseudo': '031US'},
                    'Se': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Br': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.4'},
                    'Kr': {'cutoff': 100, 'dual': 8,  'pseudo': '031US'},
                    'Rb': {'cutoff': 50,  'dual': 4,  'pseudo': 'SG15'},
                    'Sr': {'cutoff': 35,  'dual': 8,  'pseudo': '100US'},
                    'Y':  {'cutoff': 35,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Zr': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Nb': {'cutoff': 35,  'dual': 8,  'pseudo': '031PAW'},
                    'Mo': {'cutoff': 35,  'dual': 4,  'pseudo': 'SG15'},
                    'Tc': {'cutoff': 30,  'dual': 4,  'pseudo': 'SG15'},
                    'Ru': {'cutoff': 40,  'dual': 4,  'pseudo': 'SG15'},
                    'Rh': {'cutoff': 45,  'dual': 8,  'pseudo': '100PAW'},
                    'Pd': {'cutoff': 55,  'dual': 8,  'pseudo': '100PAW'},
                    'Ag': {'cutoff': 35,  'dual': 8,  'pseudo': 'GBRV-1.4'},
                    'Cd': {'cutoff': 40,  'dual': 8,  'pseudo': '031US'},
                    'In': {'cutoff': 35,  'dual': 8,  'pseudo': '031US'},
                    'Sn': {'cutoff': 35,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Sb': {'cutoff': 40,  'dual': 8,  'pseudo': 'GBRV-1.4'},
                    'Te': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'I':  {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Xe': {'cutoff': 120, 'dual': 8,  'pseudo': '100US'},
                    'Cs': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Ba': {'cutoff': 40,  'dual': 4,  'pseudo': 'SG15'},
                    'Hf': {'cutoff': 35,  'dual': 8,  'pseudo': '031US'},
                    'Ta': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'W':  {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Re': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Os': {'cutoff': 35,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Ir': {'cutoff': 40,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Pt': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.4'},
                    'Au': {'cutoff': 45,  'dual': 4,  'pseudo': 'SG15'},
                    'Hg': {'cutoff': 30,  'dual': 8,  'pseudo': 'GBRV-1.2'},
                    'Tl': {'cutoff': 30,  'dual': 8,  'pseudo': '031US'},
                    'Pb': {'cutoff': 40,  'dual': 8,  'pseudo': '031PAW'},
                    'Bi': {'cutoff': 35,  'dual': 8,  'pseudo': '031PAW'},
                    'Po': {'cutoff': 45,  'dual': 8,  'pseudo': '100US'},
                    'Rn': {'cutoff': 45,  'dual': 8,  'pseudo': '100US'},
                    'La': {'cutoff': 55,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Ce': {'cutoff': 45,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Pr': {'cutoff': 50,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Nd': {'cutoff': 40,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Sm': {'cutoff': 40,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Eu': {'cutoff': 55,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Tb': {'cutoff': 40,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Dy': {'cutoff': 40,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Ho': {'cutoff': 40,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Er': {'cutoff': 40,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Tm': {'cutoff': 40,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Yb': {'cutoff': 40,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                    'Lu': {'cutoff': 45,  'dual': 8,  'pseudo': 'Wentzcovitch'},
                }
            }

    def setup_kpoints(self):
        """
        Define the k-point mesh for the relax and scf calculations. Also get the k-point path for
        the bands calculation for the initial input structure from SeeKpath
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
            except KeyError as exception:
                self.abort_nowait('failed to retrieve the cutoff or dual factor for {}'.format(kind))

        natoms = len(structure.sites)
        conv_thr = self.ctx.protocol['convergence_threshold'] * natoms

        self.ctx.inputs['parameters'] = {
            'CONTROL': {
                'restart_mode': 'from_scratch',
                'tstress': self.ctx.protocol['tstress'],
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
        inputs = dict(self.ctx.inputs)

        options = inputs['options']
        settings = inputs['settings']
        parameters = inputs['parameters']

        # Final input preparation, wrapping dictionaries in ParameterData nodes
        inputs['kpoints_mesh'] = self.ctx.kpoints_mesh
        inputs['parameters'] = ParameterData(dict=parameters)
        inputs['settings'] = ParameterData(dict=settings)
        inputs['options'] = ParameterData(dict=options)
        inputs['relax'] = {
            'kpoints_distance': Float(self.ctx.protocol['kpoints_mesh_density']),
            'parameters': ParameterData(dict=parameters),
            'settings': ParameterData(dict=settings),
            'options': ParameterData(dict=options),
            'meta_convergence': Bool(False),
            'relaxation_scheme': Str('vc-relax'),
            'volume_convergence': Float(0.01)
        }

        running = submit(PwBandsWorkChain, **inputs)

        self.report('launching PwBandsWorkChain<{}>'.format(running.pid))

        return ToContext(workchain_bands=running)

    def run_results(self):
        """
        Attach the relevant output nodes from the band calculation to the workchain outputs
        for convenience
        """
        self.report('workchain succesfully completed')

        for link_label in ['primitive_structure', 'seekpath_parameters', 'scf_parameters', 'band_parameters', 'band_structure']:
            if link_label in self.ctx.workchain_bands.out:
                node = self.ctx.workchain_bands.get_outputs_dict()[link_label]
                self.out(link_label, node)
                self.report("attaching {}<{}> as an output node with label '{}'"
                    .format(node.__class__.__name__, node.pk, link_label))
