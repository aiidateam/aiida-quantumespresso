# -*- coding: utf-8 -*-
"""Workchain to compute all X-ray absorption spectra for a given structure.

Uses QuantumESPRESSO pw.x and xspectra.x.
"""
import warnings

from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, if_
from aiida.orm import UpfData as aiida_core_upf
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
from aiida_pseudo.data.pseudo import UpfData as aiida_pseudo_upf

from aiida_quantumespresso.utils.hubbard import HubbardStructureData
from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.workflows.protocols.utils import ProtocolMixin, recursive_merge

PwCalculation = CalculationFactory('quantumespresso.pw')
PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
PwRelaxWorkChain = WorkflowFactory('quantumespresso.pw.relax')
XyData = DataFactory('core.array.xy')
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from aiida_quantumespresso.calculations.functions.xspectra.get_spectra_by_element import get_spectra_by_element
    XspectraBaseWorkChain = WorkflowFactory('quantumespresso.xspectra.base')
    XspectraCoreWorkChain = WorkflowFactory('quantumespresso.xspectra.core')

warnings.warn(
    'This module is deprecated and will be removed soon as part of migrating XAS and XPS workflows to a new repository.'
    '\nThe new repository can be found at: https://github.com/aiidaplugins/aiida-qe-xspec.', FutureWarning
)


class XspectraCrystalWorkChain(ProtocolMixin, WorkChain):
    """Workchain to compute all X-ray absorption spectra for a given structure using Quantum ESPRESSO.

    The WorkChain follows the process required to compute all the K-edge XAS spectra for each
    element in a given structure. The WorkChain itself firstly calls the PwRelaxWorkChain to
    relax the input structure, then determines the input settings for each XAS
    calculation automatically using ``get_xspectra_structures()``:

        - Firstly the input structure is converted to its conventional standard cell using
          ``spglib`` and detects the space group number for the conventional cell.
        - Symmetry analysis of the standardized structure using ``spglib`` is then used to
          determine the number of non-equivalent atomic sites in the structure for each
          element considered for analysis.

    Using the symmetry data returned from ``get_xspectra_structures``, input structures for
    the XspectraCoreWorkChain are generated from the standardized structure by converting each
    to a supercell with cell dimensions of at least 8.0 angstroms in each periodic dimension -
    required in order to sufficiently reduce the unphysical interaction of the core-hole with
    neighbouring images. The size of the minimum size requirement can be overriden by the
    user if required. The WorkChain then uses the space group number to set the list of
    polarisation vectors for the ``XspectraCoreWorkChain`` to compute for all subsequent
    calculations.
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
            XspectraCoreWorkChain,
            namespace='core',
            exclude=(
                'kpoints', 'core_hole_pseudos', 'eps_vectors', 'structure', 'xs_plot'
            ),
            namespace_options={
                'help': ('Input parameters for the basic xspectra workflow (core-hole SCF + XAS.'),
                'validator': None
            }
        )
        spec.input_namespace(
            'core_hole_pseudos',
            # Accept both types of UpfData node
            valid_type=(aiida_core_upf, aiida_pseudo_upf),
            dynamic=True,
            help=(
                'Dynamic namespace for pairs of excited-state pseudopotentials for each absorbing'
                ' element. Must use the mapping "{element}" : {Upf}".'
            )
        )
        spec.input_namespace(
            'gipaw_pseudos',
            # Accept both types of UpfData node
            valid_type=(aiida_core_upf, aiida_pseudo_upf),
            dynamic=True,
            help=(
                'Dynamic namespace for pairs of ground-state pseudopotentials for each absorbing'
                ' element. Must use the mapping "{element}" : {Upf}.'
            )
        )
        spec.input(
            'core_hole_treatments',
            valid_type=orm.Dict,
            required=False,
            help=('Optional dictionary to set core-hole treatment to given elements present. '
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
            'elements_list',
            valid_type=orm.List,
            help=(
            'The list of elements to be considered for analysis, each must be a valid element of the periodic table.'
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
        spec.input(
            'upf2plotcore_code',
            valid_type=orm.AbstractCode,
            required=False,
            help=(
                'Code node for the upf2plotcore.sh ShellJob code.'
            )
        )
        spec.input(
            'clean_workdir',
            valid_type=orm.Bool,
            default=lambda: orm.Bool(False),
            help=('If `True`, work directories of all called calculations will be cleaned at the end of execution.'),
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
            'return_all_powder_spectra',
            valid_type=orm.Bool,
            default=lambda: orm.Bool(False),
            help=('If ``True``, the WorkChain will return all ``powder_spectrum`` nodes from each '
                  '``XspectraCoreWorkChain`` sub-process.')
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
        spec.input_namespace(
            'core_wfc_data',
            valid_type=orm.SinglefileData,
            dynamic=True,
            required=False,
            help=('Input namespace to provide core wavefunction inputs for each element. Must follow the format: '
                   '``core_wfc_data__{symbol} = {node}``')
        )
        spec.input_namespace(
            'symmetry_data',
            valid_type=(orm.Dict, orm.Int),
            dynamic=True,
            required=False,
            help=(
                'Input namespace to define equivalent sites and spacegroup number for the system. If defined, will '
                'skip symmetry analysis and structure standardization. Use *only* if symmetry data are known '
                'for certain. Requires ``spacegroup_number`` (Int) and ``equivalent_sites_data`` (Dict) to be '
                'defined separately. All keys in `equivalent_sites_data` must be formatted as "site_<site_index>". '
                'See docstring of `get_xspectra_structures` for more information about inputs.'
            )
        )
        spec.inputs.validator = cls.validate_inputs
        spec.outline(
            cls.setup,
            if_(cls.should_run_relax)(
                cls.run_relax,
                cls.inspect_relax,
            ),
            cls.get_xspectra_structures,
            if_(cls.should_run_upf2plotcore)(
                cls.run_upf2plotcore,
                cls.inspect_upf2plotcore,
            ),
            cls.run_all_xspectra_core,
            cls.inspect_all_xspectra_core,
            cls.results,
        )

        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_RELAX', message='The Relax sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_XSPECTRA', message='One or more XSpectra workflows failed')
        spec.exit_code(403, 'ERROR_NO_GIPAW_INFO_FOUND', message='The pseudos for one or more absorbing elements'
                       ' contain no GIPAW information.')
        spec.output(
            'optimized_structure',
            valid_type=orm.StructureData,
            required=False,
            help='The optimized structure from the ``relax`` process.',
        )
        spec.output(
            'standardized_structure',
            valid_type=orm.StructureData,
            required=False,
            help='The standardized crystal structure used to generate structures for XSpectra sub-processes.',
        )
        spec.output(
            'supercell_structure',
            valid_type=orm.StructureData,
            help='The supercell of ``outputs.standardized_structure`` used to generate structures for'
            ' XSpectra sub-processes.'
        )
        spec.output(
            'symmetry_analysis_data',
            valid_type=orm.Dict,
            help='The output parameters from ``get_xspectra_structures()``.'
        )
        spec.output(
            'parameters_relax',
            valid_type=orm.Dict,
            required=False,
            help='The output_parameters of the relax step.'
        )
        spec.output_namespace(
            'parameters_scf',
            valid_type=orm.Dict,
            required=False,
            dynamic=True,
            help='The output parameters of each ``PwBaseWorkChain`` performed in each ``XspectraCoreWorkChain``.'
        )
        spec.output_namespace(
            'parameters_xspectra',
            valid_type=orm.Dict,
            required=False,
            dynamic=True,
            help='The output dictionaries of each `XspectraCalculation` performed',
        )
        spec.output_namespace(
            'powder_spectra',
            valid_type=orm.XyData,
            required=False,
            dynamic=True,
            help='All the spectra generated by the WorkChain.'
        )
        spec.output_namespace(
            'final_spectra',
            valid_type=orm.XyData,
            dynamic=True,
            help='The fully-resolved spectra for each element'
        )
        # yapf: disable

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        from ..protocols import xspectra as protocols
        return files(protocols) / 'crystal.yaml'

    @classmethod
    def get_builder_from_protocol( # pylint: disable=too-many-statements
        cls, pw_code, xs_code, structure, pseudos, upf2plotcore_code=None, core_wfc_data=None,
        core_hole_treatments=None, protocol=None, overrides=None, elements_list=None,
        options=None, **kwargs
    ): # pylint: enable:too-many-statments
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param pw_code: the ``Code`` instance configured for the ``quantumespresso.pw``
            plugin.
        :param xs_code: the ``Code`` instance configured for the
            ``quantumespresso.xspectra`` plugin.
        :param upf2plotcore_code: the AiiDA-Shell ``Code`` instance configured for the
                                  upf2plotcore shell script.
        :param structure: the ``StructureData`` instance to use.
        :param pseudos: the core-hole pseudopotential pairs (ground-state and
                        excited-state) for the elements to be calculated. These must
                        use the mapping of {"element" : {"core_hole" : <upf>,"gipaw" : <upf>}}
        :param protocol: the protocol to use. If not specified, the default will be used.
        :param overrides: optional dictionary of inputs to override the defaults of the
                          XspectraWorkChain itself.
        :param kwargs: additional keyword arguments that will be passed to the
            ``get_builder_from_protocol`` of all the sub processes that are called by this
            workchain.
        :return: a process builder instance with all inputs defined ready for launch.
        """

        inputs = cls.get_protocol_inputs(protocol, overrides)

        pw_args = (pw_code, structure, protocol)

        relax = PwRelaxWorkChain.get_builder_from_protocol(
            *pw_args, overrides=inputs.get('relax', None), options=options, **kwargs
        )
        core_scf = PwBaseWorkChain.get_builder_from_protocol(
            *pw_args, overrides=inputs.get('core', {}).get('scf'), options=options, **kwargs
        )
        core_xspectra = XspectraBaseWorkChain.get_protocol_inputs(
            protocol,
            overrides=inputs.get('core', {}).get('xs_prod')
        )

        if options:
            core_xspectra['xspectra']['metadata']['options'] = recursive_merge(
                core_xspectra['xspectra']['metadata']['options'],
                options
            )

        relax.pop('clean_workdir', None)
        relax.pop('structure', None)
        relax.pop('base_final_scf', None)

        abs_atom_marker = orm.Str(inputs['abs_atom_marker'])
        # pylint: disable=no-member
        builder = cls.get_builder()
        builder.relax = relax
        builder.structure = structure
        builder.abs_atom_marker = abs_atom_marker
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])
        builder.return_all_powder_spectra = orm.Bool(inputs['return_all_powder_spectra'])
        builder.core.scf = core_scf
        builder.core.xs_prod.xspectra.code = xs_code
        builder.core.xs_prod.xspectra.parameters = orm.Dict(core_xspectra['xspectra']['parameters'])
        builder.core.xs_prod.xspectra.metadata = core_xspectra['xspectra'].get('metadata')
        builder.core.xs_prod.kpoints_distance = orm.Float(core_xspectra['kpoints_distance'])
        builder.core.get_powder_spectrum = orm.Bool(True)
        builder.core.abs_atom_marker = abs_atom_marker
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
                for element in pseudos:
                    core_hole_pseudos[element] = pseudos[element]['core_hole']
                    gipaw_pseudos[element] = pseudos[element]['gipaw']
        # if no elements list is given, we instead initalise the pseudos dict with all
        # elements in the structure. Since we require a pseudo pair for each element to
        # calculate, we generate a list of elements based on which pseudos are provided.
        else:
            builder.elements_list = orm.List(list(pseudos.keys()))
            for element in pseudos:
                core_hole_pseudos[element] = pseudos[element]['core_hole']
                gipaw_pseudos[element] = pseudos[element]['gipaw']
        builder.core_hole_pseudos = core_hole_pseudos
        builder.gipaw_pseudos = gipaw_pseudos
        if core_hole_treatments:
            builder.core_hole_treatments = orm.Dict(dict=core_hole_treatments)
        if core_wfc_data:
            builder.core_wfc_data = core_wfc_data
        elif upf2plotcore_code:
            builder.upf2plotcore_code = upf2plotcore_code
        else:
            raise ValueError(
                'No code node for upf2plotcore.sh or core wavefunction data were provided.'
            )
        # pylint: enable=no-member
        return builder


    @staticmethod
    def validate_inputs(inputs, _): # pylint: disable=too-many-return-statements
        """Validate the inputs before launching the WorkChain."""
        structure = inputs['structure']
        kinds_present = [kind.name for kind in structure.kinds]
        elements_present = sorted([kind.symbol for kind in structure.kinds])

        absorbing_elements_list = sorted(inputs['elements_list'])
        extra_elements = []
        for element in absorbing_elements_list:
            if element not in elements_present:
                extra_elements.append(element)
        if len(extra_elements) > 0:
            return (
                f'Some elements in ``elements_list`` {extra_elements} do not exist in the'
                f' structure provided {elements_present}.'
            )

        abs_atom_marker = inputs['abs_atom_marker'].value
        if abs_atom_marker in kinds_present:
            return (
                f'The marker given for the absorbing atom ("{abs_atom_marker}") matches an existing Kind in the '
                f'input structure ({kinds_present}).'
            )

        if not inputs['core']['get_powder_spectrum'].value:
            return (
                'The ``get_powder_spectrum`` input for the XspectraCoreWorkChain namespace must be ``True``.'
            )

        if 'upf2plotcore_code' not in inputs and 'core_wfc_data' not in inputs:
            return (
                'Neither a ``Code`` node for upf2plotcore.sh or a set of ``core_wfc_data`` were provided.'
            )

        if 'core_wfc_data' in inputs:
            core_wfc_data_list = sorted(inputs['core_wfc_data'].keys())
            if core_wfc_data_list != absorbing_elements_list:
                return (
                    f'The ``core_wfc_data`` provided ({core_wfc_data_list}) does not match the list of'
                    f' absorbing elements ({absorbing_elements_list})'
                )
            empty_core_wfc_data = []
            for key, value in inputs['core_wfc_data'].items():
                header_line = value.get_content()[:40]
                try:
                    num_core_states = int(header_line.split(' ')[5])
                except: # pylint: disable=bare-except
                    return (
                        'The core wavefunction data file is not of the correct format'
                    ) # pylint: enable=bare-except
                if num_core_states == 0:
                    empty_core_wfc_data.append(key)
            if len(empty_core_wfc_data) > 0:
                return (
                    f'The ``core_wfc_data`` provided for elements {empty_core_wfc_data} do not contain '
                    'any wavefunction data.'
                )

        if 'symmetry_data' in inputs:
            spacegroup_number = inputs['symmetry_data']['spacegroup_number'].value
            equivalent_sites_data = inputs['symmetry_data']['equivalent_sites_data'].get_dict()
            if spacegroup_number <= 0 or spacegroup_number >= 231:
                return (
                    f'Input spacegroup number ({spacegroup_number}) outside of valid range (1-230).'
                )

            input_elements = []
            required_keys = sorted(['symbol', 'multiplicity', 'kind_name', 'site_index'])
            invalid_entries = []
            # We check three things here: (1) are there any site indices which are outside of the possible
            # range of site indices (2) do we have all the required keys for each entry,
            # and (3) is there a mismatch between `absorbing_elements_list` and the elements specified
            # in the entries of `equivalent_sites_data`. These checks are intended only to avoid a crash.
            # We assume otherwise that the user knows what they're doing and has set everything else
            # to their preferences correctly.
            for site_label, value in equivalent_sites_data.items():
                if not set(required_keys).issubset(set(value.keys())) :
                    invalid_entries.append(site_label)
                elif value['symbol'] not in input_elements:
                    input_elements.append(value['symbol'])
                    if value['site_index'] < 0 or value['site_index'] >= len(structure.sites):
                        return (
                            f'The site index for {site_label} ({value["site_index"]}) is outside the range of '
                            + f'sites within the structure (0-{len(structure.sites) -1}).'
                        )

            if len(invalid_entries) != 0:
                return (
                    f'The required keys ({required_keys}) were not found in the following entries: {invalid_entries}'
                )

            sorted_input_elements = sorted(input_elements)
            if sorted_input_elements != absorbing_elements_list:
                return (f'Elements defined for sites in `equivalent_sites_data` ({sorted_input_elements}) '
                 f'do not match the list of absorbing elements ({absorbing_elements_list})')


    # pylint: enable=too-many-return-statements
    def setup(self):
        """Set required context variables."""
        if 'core_wfc_data' in self.inputs.keys():
            self.ctx.core_wfc_data = self.inputs.core_wfc_data


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

        optimized_structure = workchain.outputs.output_structure
        self.ctx.optimized_structure = optimized_structure
        self.out('optimized_structure', optimized_structure)

    def get_xspectra_structures(self):
        """Perform symmetry analysis of the relaxed structure and get all marked structures for XSpectra."""
        from aiida_quantumespresso.workflows.functions.get_xspectra_structures import get_xspectra_structures

        elements_list = self.inputs.elements_list

        inputs = {
            'absorbing_elements_list' : elements_list,
            'absorbing_atom_marker' : self.inputs.abs_atom_marker,
            'metadata' : {
                'call_link_label' : 'get_xspectra_structures'
            }
        }
        if 'structure_preparation_settings' in self.inputs:
            optional_cell_prep = self.inputs.structure_preparation_settings
            for key, node in optional_cell_prep.items():
                inputs[key] = node

        if isinstance(self.inputs.structure, HubbardStructureData):
            # This must be False in the case of HubbardStructureData, otherwise get_xspectra_structures will except
            inputs['standardize_structure'] = orm.Bool(False)

        if 'spglib_settings' in self.inputs:
            inputs['spglib_settings'] = self.inputs.spglib_settings

        if 'symmetry_data' in self.inputs:
            inputs['parse_symmetry'] = orm.Bool(False)
            input_sym_data = self.inputs.symmetry_data
            inputs['equivalent_sites_data'] = input_sym_data['equivalent_sites_data']
            inputs['spacegroup_number'] = input_sym_data['spacegroup_number']

        if 'relax' in self.inputs:
            result = get_xspectra_structures(self.ctx.optimized_structure, **inputs)
        else:
            result = get_xspectra_structures(self.inputs.structure, **inputs)

        supercell = result.pop('supercell')
        out_params = result.pop('output_parameters')
        spacegroup_number = out_params['spacegroup_number']
        if out_params.get_dict()['structure_is_standardized']:
            standardized = result.pop('standardized_structure')
            self.out('standardized_structure', standardized)
        if spacegroup_number in range(1, 75): # trichoric system
            self.ctx.eps_vectors = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
        if spacegroup_number in range(75, 195): # dichoric system
            self.ctx.eps_vectors = [[1., 0., 0.], [0., 0., 1.]]
        if spacegroup_number in range(195, 231): # isochoric system
            self.ctx.eps_vectors = [[1., 0., 0.]]

        structures_to_process = {f'{Key.split("_")[0]}_{Key.split("_")[1]}' : Value for Key, Value in result.items()}
        self.ctx.structures_to_process = structures_to_process
        self.ctx.equivalent_sites_data = out_params['equivalent_sites_data']

        self.out('supercell_structure', supercell)
        self.out('symmetry_analysis_data', out_params)

    def should_run_upf2plotcore(self):
        """If core wavefunction data files are specified, we skip the upf2plotcore step."""
        return 'core_wfc_data' not in self.inputs

    def run_upf2plotcore(self):
        """Run the upf2plotcore.sh utility script for each element and return the core-wavefunction data."""

        ShellJob = CalculationFactory('core.shell') # pylint: disable=invalid-name
        elements_list = self.inputs.elements_list.get_list()

        shelljobs = {}
        for element in elements_list:
            upf = self.inputs.gipaw_pseudos[f'{element}']

            shell_inputs = {}

            shell_inputs['code'] = self.inputs.upf2plotcore_code
            shell_inputs['nodes'] = {'upf': upf}
            shell_inputs['metadata'] = {
                'call_link_label': f'upf2plotcore_{element}',
                'options' : {
                    'filename_stdin' : upf.filename,
                    'resources' : {
                        'num_machines' : 1,
                        'num_mpiprocs_per_machine' : 1
                    }
                }
            }

            future_shelljob = self.submit(ShellJob, **shell_inputs)
            self.report(f'Launching upf2plotcore.sh for {element}<{future_shelljob.pk}>')
            shelljobs[f'upf2plotcore_{element}'] = future_shelljob

        return ToContext(**shelljobs)


    def inspect_upf2plotcore(self):
        """Check that the outputs from the upf2plotcore step have yielded meaningful results.

        This will simply check that the core wavefunction data returned contains at least
        one core state and return an error if this is not the case.
        """

        labels = self.inputs.elements_list.get_list()
        for label in labels:
            shelljob_node = self.ctx[f'upf2plotcore_{label}']
            core_wfc_data = shelljob_node.outputs.stdout
            header_line = core_wfc_data.get_content()[:40]
            num_core_states = int(header_line.split(' ')[5])
            if num_core_states == 0:
                return self.exit_codes.ERROR_NO_GIPAW_INFO_FOUND

    def run_all_xspectra_core(self): # pylint: disable=too-many-statements
        """Call all XspectraCoreWorkChains required to compute all requested spectra."""

        structures_to_process = self.ctx.structures_to_process
        equivalent_sites_data = self.ctx.equivalent_sites_data
        abs_atom_marker = self.inputs.abs_atom_marker.value

        xspectra_core_workchains = {}
        for site in structures_to_process:
            inputs = AttributeDict(self.exposed_inputs(XspectraCoreWorkChain, namespace='core'))
            structure = structures_to_process[site]
            inputs.structure = structure
            abs_element = equivalent_sites_data[site]['symbol']
            abs_atom_kind = equivalent_sites_data[site]['kind_name']

            if 'core_hole_treatments' in self.inputs:
                ch_treatments = self.inputs.core_hole_treatments.get_dict()
                ch_treatment = ch_treatments.get(abs_element, 'full')
            else:
                ch_treatment = 'full'

            inputs.metadata.call_link_label = f'{site}_xspectra'
            inputs.eps_vectors = orm.List(list=self.ctx.eps_vectors)
            if 'core_wfc_data' in self.inputs:
                inputs.core_wfc_data = self.inputs.core_wfc_data[abs_element]
            else:
                inputs.core_wfc_data = self.ctx[f'upf2plotcore_{abs_element}'].outputs.stdout

            # Get the given settings for the SCF inputs and then overwrite them with the
            # chosen core-hole approximation, then apply the correct pseudopotential pair.
            scf_inputs = inputs.scf.pw
            scf_params = scf_inputs.parameters.get_dict()
            ch_inputs = XspectraCoreWorkChain.get_treatment_inputs(treatment=ch_treatment)
            new_scf_params = recursive_merge(left=ch_inputs, right=scf_params)

            # Set the absorbing species index (`xiabs`) for the xspectra.x input.
            new_xs_params = inputs.xs_prod.xspectra.parameters.get_dict()
            kinds_present = sorted([kind.name for kind in structure.kinds])
            abs_species_index = kinds_present.index(abs_atom_marker) + 1
            new_xs_params['INPUT_XSPECTRA']['xiabs'] = abs_species_index

            # Set `starting_magnetization` if we are using an XCH approximation, using
            # the absorbing species as a reasonable place for the unpaired electron.
            # Alternatively, ensure the starting magnetic moment is a reasonable guess
            # given the input parameters. (e.g. it conforms to an existing magnetic
            # structure already defined for the system)

            # TODO: we need to re-visit the core-hole treatment settings,
            # in order to avoid the need for fudges like these and set these at
            # submission rather than inside the WorkChain itself.
            if 'starting_magnetization' in new_scf_params['SYSTEM']:
                inherited_mag =  new_scf_params['SYSTEM']['starting_magnetization'][abs_atom_kind]
                if ch_treatment not in ['xch_smear', 'xch_fixed']:
                    new_scf_params['SYSTEM']['starting_magnetization'][abs_atom_marker] = inherited_mag
                else: # if there is meant to be an unpaired electron, give it to the absorbing atom.
                    if inherited_mag == 0: # set it to 1, if it would be neutral in the ground-state.
                        new_scf_params['SYSTEM']['starting_magnetization'][abs_atom_marker] =  1
                    else: # assume that it takes the same magnetic configuration as the kind that it replaces.
                        new_scf_params['SYSTEM']['starting_magnetization'][abs_atom_marker] =  inherited_mag
            elif ch_treatment in ['xch_smear', 'xch_fixed']:
                new_scf_params['SYSTEM']['starting_magnetization'] = {abs_atom_marker : 1}

            # remove any duplicates created from the "core_hole_treatments.yaml" defaults
            for key in new_scf_params['SYSTEM'].keys():
                if 'starting_magnetization(' in key:
                    new_scf_params['SYSTEM'].pop(key, None)

            core_hole_pseudo = self.inputs.core_hole_pseudos[abs_element]
            gipaw_pseudo = self.inputs.gipaw_pseudos[abs_element]
            inputs.scf.pw.pseudos[abs_atom_marker] = core_hole_pseudo
            # Check how many instances of the absorbing element are present and assign
            # each the GIPAW pseudo if they are not the absorbing atom itself.
            abs_element_kinds = []
            for kind in structure.kinds:
                if kind.symbol == abs_element and kind.name != abs_atom_marker:
                    abs_element_kinds.append(kind.name)
            if len(abs_element_kinds) > 0:
                for kind_name in abs_element_kinds:
                    scf_inputs['pseudos'][kind_name] = gipaw_pseudo
            else: # if there is only one atom of the absorbing element, pop the GIPAW pseudo to avoid a crash
                scf_inputs['pseudos'].pop(abs_element, None)

            scf_inputs.parameters = orm.Dict(new_scf_params)
            inputs.scf.pw = scf_inputs
            inputs.xs_prod.xspectra.parameters = orm.Dict(new_xs_params)

            inputs = prepare_process_inputs(XspectraCoreWorkChain, inputs)

            future = self.submit(XspectraCoreWorkChain, **inputs)
            xspectra_core_workchains[site] = future
            self.report(f'launched XspectraCoreWorkChain for {site}<{future.pk}>')

        return ToContext(**xspectra_core_workchains) # pylint: enable=too-many-statements

    def inspect_all_xspectra_core(self):
        """Check that all the XspectraCoreWorkChain sub-processes finished sucessfully."""

        labels = self.ctx.structures_to_process.keys()
        work_chain_nodes = {
            label : self.ctx[label] for label in labels
        }
        failed_work_chains = []
        for label, work_chain in work_chain_nodes.items():
            if not work_chain.is_finished_ok:
                failed_work_chains.append(work_chain)
                self.report(f'XspectraCoreWorkChain for ({label}) failed with exit status {work_chain.exit_status}')
        if len(failed_work_chains) > 0:
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_XSPECTRA

    def results(self):
        """Compile all output spectra, organise and post-process all computed spectra, and send to outputs."""

        labels = self.ctx.structures_to_process.keys()
        spectra_nodes = {
            label : self.ctx[label].outputs.powder_spectrum for label in labels
        }
        spectra_nodes['metadata'] = {'call_link_label' : 'compile_final_spectra'}

        equivalent_sites_data = self.ctx.equivalent_sites_data
        elements_list = self.inputs.elements_list
        final_spectra = get_spectra_by_element(elements_list, equivalent_sites_data, **spectra_nodes)

        self.out('final_spectra', final_spectra)

        if self.inputs.return_all_powder_spectra.value:
            spectra_nodes.pop('metadata', None)
            self.out('powder_spectra', spectra_nodes)

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
