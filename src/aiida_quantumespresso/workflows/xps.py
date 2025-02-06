# -*- coding: utf-8 -*-
"""Workchain to compute the X-ray photoelectron spectroscopy (XPS) for a given structure.

Uses QuantumESPRESSO pw.x.
"""
import pathlib
from typing import Optional, Union
import warnings

from aiida import orm
from aiida.common import AttributeDict, ValidationError
from aiida.engine import ToContext, WorkChain, if_
from aiida.orm.nodes.data.base import to_aiida_type
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida_pseudo.data.pseudo import UpfData
import yaml

from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.workflows.protocols.utils import ProtocolMixin, recursive_merge

PwCalculation = CalculationFactory('quantumespresso.pw')
PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
PwRelaxWorkChain = WorkflowFactory('quantumespresso.pw.relax')
XyData = DataFactory('core.array.xy')

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from aiida_quantumespresso.calculations.functions.xspectra.get_xps_spectra import get_spectra_by_element

warnings.warn(
    'This module is deprecated and will be removed soon as part of migrating XAS and XPS workflows to a new repository.'
    '\nThe new repository can be found at: https://github.com/aiidaplugins/aiida-qe-xspec.', FutureWarning
)


def validate_inputs(inputs, _):
    """Validate the inputs before launching the WorkChain."""
    structure = inputs['structure']
    elements_present = [kind.name for kind in structure.kinds]
    abs_atom_marker = inputs['abs_atom_marker'].value
    if abs_atom_marker in elements_present:
        raise ValidationError(
            f'The marker given for the absorbing atom ("{abs_atom_marker}") matches an existing Kind in the '
            f'input structure ({elements_present}).'
        )
    if 'elements_list' in inputs:
        absorbing_elements_list = sorted(inputs['elements_list'])
        if inputs['calc_binding_energy'].value:
            ce_list = sorted(inputs['correction_energies'].get_dict().keys())
            if ce_list != absorbing_elements_list:
                raise ValidationError(
                    f'The ``correction_energies`` provided ({ce_list}) does not match the list of'
                    f' absorbing elements ({absorbing_elements_list})'
                )


class XpsWorkChain(ProtocolMixin, WorkChain):
    """Workchain to compute X-ray photoelectron spectra (XPS) for a given structure.

    The WorkChain itself firstly calls the PwRelaxWorkChain to relax the input structure if
    required. Then determines the input settings for each XPS calculation automatically using
    ``get_xspectra_structures()``. The input structures are generated from the standardized
    structure by converting each to a supercell with cell dimensions of at least 8.0 angstrom
    in each periodic dimension in order to sufficiently reduce the unphysical interaction
    of the core-hole with neighbouring images. The size of the minimum size requirement can be
    overriden by the user if required. Then the standard Delta-Self-Consistent-Field (ΔSCF)
    method is used to get the XPS binding energy. Finally, the XPS spectrum is calculated
    using the Voigt profile.

    """

    @classmethod
    def define(cls, spec):
        """Define the process specification."""

        super().define(spec)
        # yapf: disable
        spec.expose_inputs(
            PwRelaxWorkChain,
            namespace='relax',
            exclude=('structure', 'clean_workdir', 'base_final_scf'),
            namespace_options={
                'help': (
                    'Input parameters for the relax process. If not specified at all, the relaxation step is skipped.'
                ),
                'required' : False,
                'populate_defaults' : False,
            }
        )
        spec.expose_inputs(
            PwBaseWorkChain,
            namespace='ch_scf',
            exclude=('pw.structure', ),
            namespace_options={
                'help': ('Input parameters for the basic xps workflow (core-hole SCF).'),
                'validator': None
            }
        )
        spec.input_namespace(
            'core_hole_pseudos',
            valid_type=(orm.UpfData, UpfData),
            dynamic=True,
            help=(
                'Dynamic namespace for pairs of excited-state pseudopotentials for each absorbing'
                ' element. Must use the mapping "{element}" : {Upf}".'
            )
        )
        spec.input_namespace(
            'gipaw_pseudos',
            valid_type=(orm.UpfData, UpfData),
            dynamic=True,
            help=(
                'Dynamic namespace for pairs of ground-state pseudopotentials for each absorbing'
                ' element. Must use the mapping "{element}" : {Upf}".'
            )
        )
        spec.input(
            'core_hole_treatments',
            valid_type=orm.Dict,
            required=False,
            help=('Optional dictionary to set core-hole treatment to all elements present. '
                  'The default full-core-hole treatment will be used if not specified.'
                 )
        )
        spec.input(
            'structure',
            valid_type=orm.StructureData,
            help=(
                'Structure to be used for calculation.'
            )
        )
        spec.input(
            'voight_gamma',
            valid_type=orm.Float,
            default=lambda: orm.Float(0.3),
            help=(
                'The gamma parameter for the Lorenzian broadening in the Voight method.'
            )
        )
        spec.input(
            'voight_sigma',
            valid_type=orm.Float,
            default=lambda: orm.Float(0.3),
            help=(
                'The sigma parameter for the gaussian broadening in the Voight method.'
            )
        )
        spec.input(
            'abs_atom_marker',
            valid_type=orm.Str,
            default=lambda: orm.Str('X'),
            help=(
                'The name for the Kind representing the absorbing atom in the structure. '
                'Will be used in all structures generated in ``get_xspectra_structures`` step.'
            ),
        )
        spec.input_namespace(
            'structure_preparation_settings',
            valid_type=(orm.Dict, orm.Float, orm.Int, orm.Bool, orm.Str),
            dynamic=True,
            required=False,
            help=(
                'Optional settings dictionary for the ``get_xspectra_structures()`` method.'
            )
        )
        spec.input(
            'spglib_settings',
            valid_type=orm.Dict,
            required=False,
            help=(
                'Optional settings dictionary for the spglib call within ``get_xspectra_structures``.'
            )
        )
        spec.input(
            'elements_list',
            valid_type=orm.List,
            required=False,
            help=(
            'The list of elements to be considered for analysis, each must be valid elements of the periodic table.'
            )
        )
        spec.input(
            'atoms_list',
            valid_type=orm.List,
            required=False,
            help=(
            'The indices of atoms to be considered for analysis.'
            )
        )
        spec.input(
            'calc_binding_energy',
            valid_type=orm.Bool,
            default=lambda: orm.Bool(False),
            help=('If `True`, run scf calculation for the supercell.'),
        )
        spec.input(
            'correction_energies',
            valid_type=orm.Dict,
            required=False,
            help=('Optional dictionary to set the correction energy to all elements present. '
                 )
        )
        spec.input(
            'clean_workdir',
            valid_type=orm.Bool,
            default=lambda: orm.Bool(False),
            help=('If `True`, work directories of all called calculations will be cleaned at the end of execution.'),
        )
        spec.input(
            'dry_run',
            valid_type=orm.Bool,
            serializer=to_aiida_type,
            required=False,
            help='Terminate workchain steps before submitting calculations (test purposes only).'
        )
        spec.inputs.validator = validate_inputs
        spec.outline(
            cls.setup,
            if_(cls.should_run_relax)(
                cls.run_relax,
                cls.inspect_relax,
            ),
            cls.prepare_structures,
            cls.run_all_scf,
            cls.inspect_all_scf,
            cls.results,
        )

        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_RELAX', message='The Relax sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_SCF', message='The SCF Pw sub processes failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_CH_SCF', message='One or more CH_SCF Pw sub processes failed')
        spec.output(
            'optimized_structure',
            valid_type=orm.StructureData,
            required=False,
            help='The optimized structure from the ``relax`` process.',
        )
        spec.output(
            'output_parameters_relax',
            valid_type=orm.Dict,
            required=False,
            help='The output_parameters of the relax step.'
        )
        spec.output(
            'standardized_structure',
            valid_type=orm.StructureData,
            required=False,
            help='The standardized crystal structure used to generate structures for XPS sub-processes.',
        )
        spec.output(
            'supercell_structure',
            valid_type=orm.StructureData,
            required=False,
            help=('The supercell of ``outputs.standardized_structure`` used to generate structures for'
            ' XPS sub-processes.')
        )
        spec.output(
            'symmetry_analysis_data',
            valid_type=orm.Dict,
            required=False,
            help='The output parameters from ``get_xspectra_structures()``.'
        )
        spec.output(
            'output_parameters_scf',
            valid_type=orm.Dict,
            required=False,
            help='The output_parameters of the scf step.'
        )
        spec.output_namespace(
            'output_parameters_ch_scf',
            valid_type=orm.Dict,
            dynamic=True,
            help='The output parameters of each ``PwBaseWorkChain`` performed``.'
        )
        spec.output_namespace(
            'chemical_shifts',
            valid_type=orm.Dict,
            dynamic=True,
            help='All the chemical shift values for each element calculated by the WorkChain.'
        )
        spec.output_namespace(
            'binding_energies',
            valid_type=orm.Dict,
            dynamic=True,
            help='All the binding energy values for each element calculated by the WorkChain.'
        )
        spec.output_namespace(
            'final_spectra_cls',
            valid_type=orm.XyData,
            dynamic=True,
            help='The fully-resolved spectra for each element based on chemical shift.'
        )
        spec.output_namespace(
            'final_spectra_be',
            valid_type=orm.XyData,
            dynamic=True,
            help='The fully-resolved spectra for each element based on binding energy.'
        )
        # yapf: disable

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        from . import protocols  # pylint: disable=relative-beyond-top-level

        # import protocols  # pylint: disable=relative-beyond-top-level
        return files(protocols) / 'xps.yaml'

    @classmethod
    def get_default_treatment(cls) -> str:
        """Return the default core-hole treatment.

        :param cls: the workflow class.
        :return: the default core-hole treatment
        """

        return cls._load_treatment_file()['default_treatment']

    @classmethod
    def get_available_treatments(cls) -> dict:
        """Return the available core-hole treatments.

        :param cls: the workflow class.
        :return: dictionary of available treatments, where each key is a treatment and value
                 is another dictionary that contains at least the key `description` and
                 optionally other keys with supplimentary information.
        """
        data = cls._load_treatment_file()
        return {treatment: {'description': values['description']} for treatment, values in data['treatments'].items()}

    @classmethod
    def get_treatment_inputs(
        cls,
        treatment: Optional[dict] = None,
        overrides: Union[dict, pathlib.Path, None] = None,
    ) -> dict:
        """Return the inputs for the given workflow class and core-hole treatment.

        :param cls: the workflow class.
        :param treatment: optional specific treatment, if not specified, the default will be used
        :param overrides: dictionary of inputs that should override those specified by the treatment. The mapping should
            maintain the exact same nesting structure as the input port namespace of the corresponding workflow class.
        :return: mapping of inputs to be used for the workflow class.
        """
        data = cls._load_treatment_file()
        treatment = treatment or data['default_treatment']

        try:
            treatment_inputs = data['treatments'][treatment]
        except KeyError as exception:
            raise ValueError(
                f'`{treatment}` is not a valid treatment. '
                'Call ``get_available_treatments`` to show available treatments.'
            ) from exception
        inputs = recursive_merge(data['default_inputs'], treatment_inputs)
        inputs.pop('description')

        if isinstance(overrides, pathlib.Path):
            with overrides.open() as file:
                overrides = yaml.safe_load(file)

        if overrides:
            return recursive_merge(inputs, overrides)

        return inputs

    @classmethod
    def _load_treatment_file(cls) -> dict:
        """Return the contents of the core-hole treatment file."""
        with cls.get_treatment_filepath().open() as file:
            return yaml.safe_load(file)

    @classmethod
    def get_treatment_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the core-hole treatments for the SCF step."""
        from importlib_resources import files

        from . import protocols

        # import protocols
        return files(protocols) / 'core_hole_treatments.yaml'

    @classmethod
    def get_builder_from_protocol(
        cls, code, structure, pseudos, core_hole_treatments=None, protocol=None,
        overrides=None, elements_list=None, atoms_list=None, options=None,
        structure_preparation_settings=None, correction_energies=None, **kwargs
    ): # pylint: disable=too-many-statements
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param structure: the ``StructureData`` instance to use.
        :param pseudos: the core-hole pseudopotential pairs (ground-state and
                        excited-state) for the elements to be calculated. These must
                        use the mapping of {"element" : {"core_hole" : <upf>, "gipaw" : <upf>}}
        :param protocol: the protocol to use. If not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the
                          XpsWorkChain itself.
        :param kwargs: additional keyword arguments that will be passed to the
            ``get_builder_from_protocol`` of all the sub processes that are called by this
            workchain.
        :return: a process builder instance with all inputs defined ready for launch.
        """

        inputs = cls.get_protocol_inputs(protocol, overrides)
        pw_args = (code, structure, protocol)
        # xspectra_args = (pw_code, xs_code, structure, protocol, upf2plotcore_code)

        relax = PwRelaxWorkChain.get_builder_from_protocol(
            *pw_args, overrides=inputs.get('relax', None), options=options, **kwargs
        )
        ch_scf = PwBaseWorkChain.get_builder_from_protocol(
            *pw_args, overrides=inputs.get('ch_scf', None), options=options, **kwargs
        )

        relax.pop('clean_workdir', None)
        relax.pop('structure', None)
        relax.pop('base_final_scf', None)
        ch_scf.pop('clean_workdir', None)
        ch_scf.pop('structure', None)

        abs_atom_marker = orm.Str(inputs['abs_atom_marker'])
        # pylint: disable=no-member
        builder = cls.get_builder()
        builder.relax = relax
        builder.ch_scf = ch_scf
        builder.structure = structure
        builder.abs_atom_marker = abs_atom_marker
        if correction_energies:
            builder.correction_energies = orm.Dict(correction_energies)
            builder.calc_binding_energy = orm.Bool(True)
        else:
            builder.calc_binding_energy = orm.Bool(False)
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])
        core_hole_pseudos = {}
        gipaw_pseudos = {}
        if elements_list:
            elements_not_present = []
            elements_present = [kind.symbol for kind in structure.kinds]
            for element in elements_list:
                if element not in elements_present:
                    elements_not_present.append(element)

            if len(elements_not_present) > 0:
                raise ValueError(
                    f'The following elements: {elements_not_present} are not present in the'
                    f' structure ({elements_present}) provided.'
                    )
            else:
                builder.elements_list = orm.List(elements_list)
                for element in elements_list:
                    core_hole_pseudos[element] = pseudos[element]['core_hole']
                    gipaw_pseudos[element] = pseudos[element]['gipaw']
        elif atoms_list:
            builder.atoms_list = orm.List(atoms_list)
            for index in atoms_list:
                element = structure.sites[index].kind_name
                core_hole_pseudos[element] = pseudos[element]['core_hole']
                gipaw_pseudos[element] = pseudos[element]['gipaw']
        # if no elements list is given, we instead initalise the pseudos dict with all
        # elements in the structure
        else:
            for element in pseudos:
                core_hole_pseudos[element] = pseudos[element]['core_hole']
                gipaw_pseudos[element] = pseudos[element]['gipaw']
        builder.core_hole_pseudos = core_hole_pseudos
        builder.gipaw_pseudos = gipaw_pseudos
        if core_hole_treatments:
            builder.core_hole_treatments = orm.Dict(dict=core_hole_treatments)
        # for get_xspectra_structures
        if structure_preparation_settings:
            builder.structure_preparation_settings = structure_preparation_settings
            if structure_preparation_settings.get('is_molecule_input').value:
                builder.ch_scf.pw.parameters.base.attributes.all['SYSTEM']['assume_isolated']='mt'
                builder.ch_scf.pw.settings=orm.Dict(dict={'gamma_only':True})
                # To ensure compatibility with the gamma_only setting, the k-points must be configured to [1, 1, 1].
                kpoints_mesh = DataFactory('core.array.kpoints')()
                kpoints_mesh.set_kpoints_mesh([1, 1, 1])
                builder.ch_scf.kpoints = kpoints_mesh
                builder.relax.base.pw.settings=orm.Dict(dict={'gamma_only':True})
        # pylint: enable=no-member
        return builder


    def setup(self):
        """Init required context variables."""
        elements_list = self.inputs.get('elements_list', None)
        atoms_list = self.inputs.get('atoms_list', None)
        if elements_list:
            self.ctx.elements_list = elements_list.get_list()
            self.ctx.atoms_list = None
        elif atoms_list:
            self.ctx.atoms_list = atoms_list.get_list()
            self.ctx.elements_list = None
        else:
            structure = self.inputs.structure
            self.ctx.elements_list = [Kind.symbol for Kind in structure.kinds]



    def should_run_relax(self):
        """If the 'relax' input namespace was specified, we relax the input structure."""
        return 'relax' in self.inputs

    def run_relax(self):
        """Run the PwRelaxWorkChain to run a relax PwCalculation."""
        inputs = AttributeDict(self.exposed_inputs(PwRelaxWorkChain, namespace='relax'))
        inputs.metadata.call_link_label = 'relax'
        inputs.structure = self.inputs.structure

        running = self.submit(PwRelaxWorkChain, **inputs)

        self.report(f'launching PwRelaxWorkChain<{running.pk}>')

        return ToContext(relax_workchain=running)

    def inspect_relax(self):
        """Verify that the PwRelaxWorkChain finished successfully."""
        workchain = self.ctx.relax_workchain

        if not workchain.is_finished_ok:
            self.report(f'PwRelaxWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX

        relaxed_structure = workchain.outputs.output_structure
        relax_params = workchain.outputs.output_parameters
        self.ctx.relaxed_structure = relaxed_structure
        self.out('optimized_structure', relaxed_structure)
        self.out('output_parameters_relax', relax_params)

    def prepare_structures(self):
        """Get a marked structure for each site.

        Analyses the given structure using ``get_xspectra_structures()`` to obtain both the
        conventional standard form of the crystal structure and the list of symmetrically
        non-equivalent sites. This list is then used to produce a supercell of the
        standardized structure for each selected site with the site marked using a unique
        label.

        If provided, this step will use inputs from ``inputs.structure_preparation_settings``
        and apply them to the CalcFunction call. The accepted inputs (format, default) are:
        - supercell_min_parameter (float, 8.0)
        - standardize_structure (bool, True)
        - is_molecule_input (bool, False)

        Input settings for the spglib analysis within ``get_xspectra_structures`` can be
        provided via ``inputs.spglib_settings`` in the form of a Dict node and must be
        formatted as {<variable_name> : <parameter>} for each variable in the
        ``get_symmetry_dataset()`` method.
        """
        from aiida_quantumespresso.workflows.functions.get_marked_structures import get_marked_structures
        from aiida_quantumespresso.workflows.functions.get_xspectra_structures import get_xspectra_structures

        input_structure = self.inputs.structure if 'relax' not in self.inputs else self.ctx.relaxed_structure
        if self.ctx.elements_list:
            elements_list = orm.List(self.ctx.elements_list)
            inputs = {
                'absorbing_elements_list' : elements_list,
                'absorbing_atom_marker' : self.inputs.abs_atom_marker,
                'metadata' : {
                    'call_link_label' : 'get_xspectra_structures'
                }
            } # populate this further once the schema for WorkChain options is figured out
            if 'structure_preparation_settings' in self.inputs:
                optional_cell_prep = self.inputs.structure_preparation_settings
                for key, node in optional_cell_prep.items():
                    inputs[key] = node
            if 'spglib_settings' in self.inputs:
                spglib_settings = self.inputs.spglib_settings
                inputs['spglib_settings'] = spglib_settings
            else:
                spglib_settings = None

            result = get_xspectra_structures(input_structure, **inputs)

            supercell = result.pop('supercell')
            out_params = result.pop('output_parameters')
            if out_params.get_dict().get('structure_is_standardized', None):
                standardized = result.pop('standardized_structure')
                self.out('standardized_structure', standardized)

            # structures_to_process = {Key : Value for Key, Value in result.items()}
            for site in ['output_parameters', 'supercell', 'standardized_structure']:
                result.pop(site, None)
            self.ctx.supercell = supercell
            self.ctx.equivalent_sites_data = out_params['equivalent_sites_data']
            self.out('supercell_structure', supercell)
            self.out('symmetry_analysis_data', out_params)
        elif self.ctx.atoms_list:
            atoms_list = orm.List(self.ctx.atoms_list)
            inputs = {
                'atoms_list' : atoms_list,
                'marker' : self.inputs.abs_atom_marker,
                'metadata' : {
                    'call_link_label' : 'get_marked_structures'
                }
            }
            result = get_marked_structures(input_structure, **inputs)
            self.ctx.supercell = input_structure
            self.ctx.equivalent_sites_data = result.pop('output_parameters').get_dict()
        structures_to_process = {f'{Key.split("_")[0]}_{Key.split("_")[1]}' : Value for Key, Value in result.items()}
        self.report(f'structures_to_process: {structures_to_process}')
        self.ctx.structures_to_process = structures_to_process

    def should_run_gs_scf(self):
        """If the 'calc_binding_energy' input namespace is True, we run a scf calculation for the supercell."""
        return self.inputs.calc_binding_energy

    def run_gs_scf(self):
        """Call ``PwBaseWorkChain`` to compute total energy for the supercell."""

        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='ch_scf'))
        inputs.pw.structure = self.ctx.supercell
        inputs.metadata.call_link_label = 'supercell_xps'

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        # pseudos for all elements to be calculated should be replaced
        for site in self.ctx.equivalent_sites_data:
            abs_element = self.ctx.equivalent_sites_data[site]['symbol']
            inputs.pw.pseudos[abs_element] = self.inputs.gipaw_pseudos[abs_element]
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report(f'launched PwBaseWorkChain for supercell<{running.pk}>')

        return running

    def inspect_scf(self):
        """Verify that the PwBaseWorkChain finished successfully."""
        workchain = self.ctx.scf_workchain

        if not workchain.is_finished_ok:
            self.report(f'PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF

        scf_params = workchain.outputs.output_parameters
        self.out('output_parameters_scf', scf_params)

    def run_all_scf(self):
        """Call all PwBaseWorkChain's required to compute total energies for each absorbing atom site."""

        # scf for supercell
        futures = {}
        if self.inputs.calc_binding_energy:
            gs_future = self.run_gs_scf()
            futures['ground_state'] = gs_future
        # scf for core hole
        structures_to_process = self.ctx.structures_to_process
        equivalent_sites_data = self.ctx.equivalent_sites_data
        abs_atom_marker = self.inputs.abs_atom_marker.value

        for site in structures_to_process:
            inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='ch_scf'))
            structure = structures_to_process[site]
            inputs.pw.structure = structure
            abs_element = equivalent_sites_data[site]['symbol']

            if 'core_hole_treatments' in self.inputs:
                ch_treatments = self.inputs.core_hole_treatments.get_dict()
                ch_treatment = ch_treatments.get(abs_element, 'xch_smear')
            else:
                ch_treatment = 'xch_smear'


            inputs.metadata.call_link_label = f'{site}_xps'

            # Get the given settings for the SCF inputs and then overwrite them with the
            # chosen core-hole approximation, then apply the correct pseudopotential pair
            scf_params = inputs.pw.parameters.get_dict()
            ch_treatment_inputs = self.get_treatment_inputs(treatment=ch_treatment)

            new_scf_params = recursive_merge(left=scf_params, right=ch_treatment_inputs)
            if ch_treatment == 'xch_smear':
                structure_kinds = [kind.name for kind in structure.kinds]
                structure_kinds.sort()
                abs_species = structure_kinds.index(abs_atom_marker)
                new_scf_params['SYSTEM'][f'starting_magnetization({abs_species + 1})'] = 1

            core_hole_pseudo = self.inputs.core_hole_pseudos[abs_element]
            inputs.pw.pseudos[abs_atom_marker] = core_hole_pseudo
            # pseudos for all elements to be calculated should be replaced
            for key in self.ctx.equivalent_sites_data:
                abs_element = self.ctx.equivalent_sites_data[key]['symbol']
                inputs.pw.pseudos[abs_element] = self.inputs.gipaw_pseudos[abs_element]
            # remove pseudo if the only element is replaced by the marker
            inputs.pw.pseudos = {kind.name: inputs.pw.pseudos[kind.name] for kind in structure.kinds}

            inputs.pw.parameters = orm.Dict(dict=new_scf_params)

            inputs = prepare_process_inputs(PwBaseWorkChain, inputs)

            future = self.submit(PwBaseWorkChain, **inputs)
            futures[site] = future
            self.report(f'launched PwBaseWorkChain for {site}<{future.pk}>')

        return ToContext(**futures)

    def inspect_all_scf(self):
        """Check that all the PwBaseWorkChain sub-processes finished sucessfully."""

        labels = self.ctx.structures_to_process.keys()
        work_chains = [self.ctx[label] for label in labels]
        failed_work_chains = []
        for work_chain, label in zip(work_chains, labels):
            if not work_chain.is_finished_ok:
                failed_work_chains.append(work_chain)
                self.report(f'PwBaseWorkChain for ({label}) failed with exit status {work_chain.exit_status}')
        if len(failed_work_chains) > 0:
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_CH_SCF

    def results(self):
        """Compile all output spectra, organise and post-process all computed spectra, and send to outputs."""

        import copy  # pylint: disable=unused-import

        site_labels = list(self.ctx.structures_to_process.keys())
        output_params_ch_scf = {label : self.ctx[label].outputs.output_parameters for label in site_labels}
        self.out('output_parameters_ch_scf', output_params_ch_scf)

        kwargs = output_params_ch_scf.copy()
        if self.inputs.calc_binding_energy:
            kwargs['ground_state'] = self.ctx['ground_state'].outputs.output_parameters
            kwargs['correction_energies'] = self.inputs.correction_energies
        kwargs['metadata'] = {'call_link_label' : 'compile_final_spectra'}

        if self.ctx.elements_list:
            elements_list = orm.List(list=self.ctx.elements_list)
        else:
            symbols = {value['symbol'] for value in self.ctx.equivalent_sites_data.values()}
            elements_list = orm.List(list(symbols))
        voight_gamma = self.inputs.voight_gamma
        voight_sigma = self.inputs.voight_sigma

        equivalent_sites_data = orm.Dict(dict=self.ctx.equivalent_sites_data)
        result = get_spectra_by_element(
            elements_list,
            equivalent_sites_data,
            voight_gamma,
            voight_sigma,
            **kwargs
        )
        final_spectra_cls = {}
        final_spectra_be = {}
        chemical_shifts = {}
        binding_energies = {}
        for key, value in result.items():
            if key.endswith('cls_spectra'):
                final_spectra_cls[key] = value
            elif key.endswith('be_spectra'):
                final_spectra_be[key] = value
            elif key.endswith('cls'):
                chemical_shifts[key] = value
            elif key.endswith('be'):
                binding_energies[key] = value
        self.out('chemical_shifts', chemical_shifts)
        self.out('final_spectra_cls', final_spectra_cls)
        if self.inputs.calc_binding_energy:
            self.out('binding_energies', binding_energies)
            self.out('final_spectra_be', final_spectra_be)


    def on_terminated(self):
        """Clean the working directories of all child calculations if ``clean_workdir=True`` in the inputs."""

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
