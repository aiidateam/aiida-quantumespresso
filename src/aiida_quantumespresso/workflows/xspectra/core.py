# -*- coding: utf-8 -*-
"""Workchain to compute the X-ray absorption spectrum for a given structure.

Uses QuantumESPRESSO pw.x and xspectra.x, requires ``aiida-shell`` to run ``upf2plotcore.sh``.
"""
import pathlib
from typing import Optional, Union
import warnings

from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, append_, if_
from aiida.orm.nodes.data.base import to_aiida_type
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory
import yaml

from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.workflows.protocols.utils import ProtocolMixin, recursive_merge

PwCalculation = CalculationFactory('quantumespresso.pw')
PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from aiida_quantumespresso.calculations.functions.xspectra.get_powder_spectrum import get_powder_spectrum
    from aiida_quantumespresso.calculations.functions.xspectra.merge_spectra import merge_spectra
    XspectraBaseWorkChain = WorkflowFactory('quantumespresso.xspectra.base')

HubbardStructureData = DataFactory('quantumespresso.hubbard_structure')

warnings.warn(
    'This module is deprecated and will be removed soon as part of migrating XAS and XPS workflows to a new repository.'
    '\nThe new repository can be found at: https://github.com/aiidaplugins/aiida-qe-xspec.', FutureWarning
)


class XspectraCoreWorkChain(ProtocolMixin, WorkChain):
    """Workchain to compute X-ray absorption spectra for a given structure using Quantum ESPRESSO.

    The workflow follows the process required to compute the XAS of an input structure: an SCF calculation is performed
    using the provided structure, which is then followed by the calculation of the XAS itself by XSpectra. The
    calculations performed by the WorkChain in a typical run will be:

    - PwSCF calculation with pw.x of the input structure with a core-hole present.
    - Generation of core-wavefunction data with upf2plotcore.sh (if requested).
    - XAS calculation with xspectra.x to compute the Lanczos coefficients and print the XANES spectra for the
      polarisation vectors requested in the input.
    - Collation of output data from pw.x and xspectra.x calculations, including a combination of XANES dipole spectra
      based on polarisation vectors to represent the powder spectrum of the structure (if requested).

    If ``run_replot = True`` is set in the inputs (defaults to False), the WorkChain will run a second xspectra.x
    calculation which replots the spectra produced from the ``xs_prod`` step. This option can be very useful for
    obtaining a final spectrum at low levels of broadening (relative to the default of 0.5 eV), particularly as higher
    levels of broadening significantly speed up the convergence of the Lanczos procedure. Inputs for the replot
    calculation are found in the ``xs_plot`` namespace.

    The core-wavefunction plot derived from the ground-state of the absorbing element can be provided as a top-level
    input or produced by the WorkChain. If left to the WorkChain, the ground-state pseudopotential assigned to the
    absorbing element will be used to generate this data using the upf2plotcore.sh utility script (via the
    ``aiida-shell`` plugin).

    In its current stage of development, the workflow requires the following:

    - An input structure where the desired absorbing atom in the system is marked as a separate Kind. The default
      behaviour for the WorkChain is to set the Kind name as 'X', however this can be changed via the `overrides`
      dictionary.
    - A code node for ``upf2plotcore``, configured for the ``aiida-shell`` plugin
      (https://github.com/sphuber/aiida-shell). Alternatively, a ``SinglefileData`` node from a previous ``ShellJob``
      run can be supplied under ``inputs.core_wfc_data``.
    - A suitable pair of pseudopotentials for the element type of the absorbing atom, one for the ground-state occupancy
      which contains GIPAW informtation for the core level of interest for the XAS (e.g. 1s in the case of a K-edge
      calculation) and the other containing a core hole. (For the moment this can be passed either via the
      ``core_hole_pseudos`` field in ``get_builder_from_protocol`` or via the overrides, but will be changed later once
      full families of core-hole pseudopotentials become available).
    """

    # pylint: disable=too-many-public-methods, too-many-statements

    @classmethod
    def define(cls, spec):
        """Define the process specification."""

        super().define(spec)
        spec.expose_inputs(
            PwBaseWorkChain,
            namespace='scf',
            exclude=('pw.parent_folder', 'pw.structure', 'clean_workdir'),
            namespace_options={
                'help': ('Input parameters for the `pw.x` calculation.'),
                'validator': cls.validate_scf,
            }
        )
        spec.expose_inputs(
            XspectraBaseWorkChain,
            namespace='xs_prod',
            exclude=('clean_workdir', 'xspectra.parent_folder', 'xspectra.core_wfc_data'),
            namespace_options={
                'help': ('Input parameters for the `xspectra.x` calculation'
                         ' to compute the Lanczos.')
            }
        )
        spec.expose_inputs(
            XspectraBaseWorkChain,
            namespace='xs_plot',
            exclude=('clean_workdir', 'xspectra.parent_folder', 'xspectra.core_wfc_data'),
            namespace_options={
                'help': ('Input parameters for the re-plot `xspectra.x` calculation of the Lanczos.'),
                'required': False,
                'populate_defaults': False
            }
        )
        spec.inputs.validator = cls.validate_inputs
        spec.input(
            'structure',
            valid_type=(orm.StructureData, HubbardStructureData),
            help=(
                'Structure to be used for calculation, with at least one site containing the `abs_atom_marker` '
                'as the kind label.'
            )
        )
        spec.input(
            'eps_vectors',
            valid_type=orm.List,
            help=(
                'The list of 3-vectors to use in XSpectra sub-processes. '
                'The number of sub-lists will subsequently define the number of XSpectra calculations to perform'
            ),
        )
        spec.input(
            'abs_atom_marker',
            valid_type=orm.Str,
            required=False,
            help=(
                'The name for the Kind representing the absorbing atom in the structure. '
                'Must corespond to a Kind within the StructureData node supplied to the calculation.'
            ),
        )
        spec.input(
            'get_powder_spectrum',
            valid_type=orm.Bool,
            default=lambda: orm.Bool(False),
            help=(
                'If `True`, the WorkChain will combine XANES dipole spectra computed using the XAS basis vectors'
                ' defined according to the `get_powder_spectrum` CalcFunction.'
            ),
        )
        spec.input(
            'core_wfc_data',
            valid_type=orm.SinglefileData,
            required=False,
            help='The core wavefunction data file extracted from the ground-state pseudo for the absorbing atom.'
        )
        spec.input(
            'run_replot',
            valid_type=orm.Bool,
            serializer=to_aiida_type,
            default=lambda: orm.Bool(False),
        )
        spec.input(
            'upf2plotcore_code',
            valid_type=orm.AbstractCode,
            required=False,
            help='The code node required for upf2plotcore.sh configured for ``aiida-shell``. '
            'Must be provided if `core_wfc_data` is not provided.'
        )
        spec.input(
            'clean_workdir',
            valid_type=orm.Bool,
            serializer=to_aiida_type,
            default=lambda: orm.Bool(False),
            help=('If `True`, work directories of all called calculation will be cleaned at the end of execution.'),
        )
        spec.input(
            'dry_run',
            valid_type=orm.Bool,
            serializer=to_aiida_type,
            required=False,
            help='Terminate workchain steps before submitting calculations (test purposes only).'
        )
        spec.outline(
            cls.setup,
            cls.run_scf,
            cls.inspect_scf,
            if_(cls.should_run_upf2plotcore)(
                cls.run_upf2plotcore,
                cls.inspect_upf2plotcore,
            ),
            cls.run_all_xspectra_prod,
            cls.inspect_all_xspectra_prod,
            if_(cls.should_run_replot)(
                cls.run_all_xspectra_plot,
                cls.inspect_all_xspectra_plot,
            ),
            cls.results,
        )

        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_SCF', message='The SCF sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_XSPECTRA', message='One or more XSpectra sub processes failed')
        spec.exit_code(
            403,
            'ERROR_NO_GIPAW_INFO_FOUND',
            message='The pseudo for the absorbing element contains no'
            ' GIPAW information.'
        )
        spec.output(
            'parameters_scf', valid_type=orm.Dict, help='The output parameters of the SCF'
            ' `PwBaseWorkChain`.'
        )
        spec.output_namespace(
            'parameters_xspectra',
            valid_type=orm.Dict,
            help='The output dictionaries of each `XspectraBaseWorkChain` performed',
            dynamic=True
        )
        spec.output(
            'spectra',
            valid_type=orm.XyData,
            help='An XyData node containing all the final spectra produced by the WorkChain.'
        )
        spec.output('powder_spectrum', valid_type=orm.XyData, required=False, help='The simulated powder spectrum')

    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        from ..protocols import xspectra as protocols
        return files(protocols) / 'core.yaml'

    @classmethod
    def get_treatment_filepath(cls) -> pathlib.Path:
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the core-hole treatments for the SCF step."""
        from importlib_resources import files

        from .. import protocols
        return files(protocols) / 'core_hole_treatments.yaml'

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
    def get_builder_from_protocol(
        cls,
        pw_code,
        xs_code,
        structure,
        upf2plotcore_code=None,
        core_wfc_data=None,
        core_hole_pseudos=None,
        core_hole_treatment=None,
        protocol=None,
        overrides=None,
        options=None,
        **kwargs
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param pw_code: the ``Code`` instance configured for the ``quantumespresso.pw``
                        plugin.
        :param xs_code: the ``Code`` instance configured for the
                        ``quantumespresso.xspectra`` plugin.
        :param structure: the ``StructureData`` instance to use.
        :param upf2plotcore_code: the aiida-shell ``Code`` instance configured for the
                                  upf2plotcore.sh shell script, used to generate the core
                                  wavefunction data.
        :param core_wfc_data:
        :param core_hole_pseudos: the core-hole pseudopotential pair (ground-state and
                                  excited-state) for the chosen absorbing atom.
        :param protocol: the protocol to use. If not specified, the default will be used.
        :param core_hole_treatment: the core-hole treatment desired for the SCF calculation,
                                    using presets found in ``core_hole_treatments.yaml``.
                                    Defaults to "full". Overrides the settings derived from
                                    the ``PwBaseWorkChain`` protocol, but is itself overriden
                                    by the ``overrides`` dictionary.
        :param overrides: optional dictionary of inputs to override the defaults of the
                          XspectraBaseWorkChain itself.
        :param options: a dictionary of options that will be recursively set for the #
                        ``metadata.options`` input of all the ``CalcJobs`` that are nested in
                        this work chain.
        :param run_replot: a bool parameter to request inputs for the re-plot step.
        :param kwargs: additional keyword arguments that will be passed to the
            ``get_builder_from_protocol`` of all the sub processes that are called by this
            workchain.
        :return: a process builder instance with all inputs defined ready for launch.
        """

        inputs = cls.get_protocol_inputs(protocol, overrides)
        pw_inputs = PwBaseWorkChain.get_protocol_inputs(protocol=protocol, overrides=inputs.get('scf', {}))
        pw_params = pw_inputs['pw']['parameters']
        kinds_present = sorted([kind.name for kind in structure.kinds])
        # Get the default inputs from the PwBaseWorkChain and override them with those
        # required for the chosen core-hole treatment
        pw_params = recursive_merge(
            left=pw_params,
            right=cls.get_treatment_inputs(
                treatment=core_hole_treatment, overrides=inputs.get('scf', {}).get('pw', {}).get('parameters', None)
            )
        )

        pw_inputs['pw']['parameters'] = pw_params
        pw_args = (pw_code, structure, protocol)

        scf = PwBaseWorkChain.get_builder_from_protocol(*pw_args, overrides=pw_inputs, options=options, **kwargs)

        scf.pop('clean_workdir', None)
        scf['pw'].pop('structure', None)

        # pylint: disable=no-member
        builder = cls.get_builder()
        builder.scf = scf

        xs_prod_inputs = XspectraBaseWorkChain.get_protocol_inputs(protocol, inputs.get('xs_prod'))
        xs_prod_parameters = xs_prod_inputs['xspectra']['parameters']
        xs_prod_metadata = xs_prod_inputs['xspectra']['metadata']
        if options:
            xs_prod_metadata['options'] = recursive_merge(xs_prod_metadata['options'], options)

        abs_atom_marker = inputs['abs_atom_marker']
        xs_prod_parameters['INPUT_XSPECTRA']['xiabs'] = kinds_present.index(abs_atom_marker) + 1
        if core_hole_pseudos:
            abs_element_kinds = []
            for kind in structure.kinds:
                if kind.name == abs_atom_marker:
                    abs_element = kind.symbol
            for kind in structure.kinds:  # run a second pass to check for multiple kinds of the same absorbing element
                if kind.symbol == abs_element and kind.name != abs_atom_marker:
                    abs_element_kinds.append(kind.name)
            builder.scf.pw.pseudos[abs_atom_marker] = core_hole_pseudos[abs_atom_marker]
            for kind_name in abs_element_kinds:
                builder.scf.pw.pseudos[kind_name] = core_hole_pseudos[abs_element]

        builder.xs_prod.xspectra.code = xs_code
        builder.xs_prod.xspectra.parameters = orm.Dict(xs_prod_parameters)
        builder.xs_prod.xspectra.metadata = xs_prod_metadata
        if xs_prod_inputs['kpoints_distance']:
            builder.xs_prod.kpoints_distance = orm.Float(xs_prod_inputs['kpoints_distance'])
        elif xs_prod_inputs['kpoints']:
            builder.xs_prod.kpoints = xs_prod_inputs['kpoints']

        if upf2plotcore_code:
            builder.upf2plotcore_code = upf2plotcore_code
        elif core_wfc_data:
            builder.core_wfc_data = core_wfc_data
        else:
            raise ValueError(
                'Either a code node for upf2plotcore.sh or an already-generated core-wavefunction'
                ' file must be given.'
            )

        builder.structure = structure
        builder.eps_vectors = orm.List(list=inputs['eps_vectors'])
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])
        builder.get_powder_spectrum = orm.Bool(inputs['get_powder_spectrum'])
        builder.abs_atom_marker = orm.Str(abs_atom_marker)
        if inputs['run_replot']:
            builder.run_replot = orm.Bool(inputs['run_replot'])
            xs_plot_inputs = XspectraBaseWorkChain.get_protocol_inputs('replot')
            xs_plot_parameters = xs_plot_inputs['xspectra']['parameters']
            xs_plot_metadata = xs_plot_inputs['xspectra']['metadata']
            if options:
                xs_plot_metadata['options'] = recursive_merge(xs_plot_metadata['options'], options)
            builder.xs_plot.xspectra.code = xs_code
            builder.xs_plot.xspectra.parameters = orm.Dict(xs_plot_parameters)
            builder.xs_plot.xspectra.metadata = xs_plot_metadata
        else:
            builder.pop('run_replot', None)
        # pylint: enable=no-member
        return builder

    @staticmethod
    def validate_scf(inputs, _):
        """Validate the scf parameters."""
        parameters = inputs['pw']['parameters'].get_dict()
        if parameters.get('CONTROL', {}).get('calculation', 'scf') != 'scf':
            return '`CONTROL.calculation` in `scf.pw.parameters` is not set to `scf`.'

    @staticmethod
    def validate_inputs(inputs, _):
        """Validate the inputs before launching the WorkChain."""

        eps_vector_list = inputs['eps_vectors'].get_list()
        if len(eps_vector_list) == 0:
            return '`eps_vectors` list empty.'
        if 'core_wfc_data' not in inputs and 'upf2plotcore_code' not in inputs:
            return 'Either a core wavefunction file or a code node for upf2plotcore.sh must be provided.'
        structure = inputs['structure']
        kinds_present = structure.kinds
        abs_atom_found = False
        for kind in kinds_present:
            if kind.name == inputs['abs_atom_marker'].value:
                abs_atom_found = True
        if not abs_atom_found:
            return (
                f'Error: the marker given for the absorbing atom ("{inputs["abs_atom_marker"].value}") ' +
                'does not appear in the structure provided.'
            )

    def setup(self):
        """Initialize context variables that are used during the logical flow of the workchain."""

        self.ctx.dry_run = 'dry_run' in self.inputs and self.inputs.dry_run.value
        self.ctx.all_lanczos_computed = False
        self.ctx.finished_lanczos = []
        self.ctx.finished_replots = []
        abs_atom_marker = self.inputs.abs_atom_marker
        structure = self.inputs.structure
        for kind in structure.kinds:
            if kind.name == abs_atom_marker:
                abs_kind = kind

        self.ctx.abs_kind = abs_kind
        if 'core_wfc_data' in self.inputs:
            self.ctx.core_wfc_data = self.inputs.core_wfc_data

    def run_scf(self):
        """Run an SCF calculation as a first step."""

        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, 'scf'))

        inputs.metadata.call_link_label = 'scf'
        inputs.pw.structure = self.inputs.structure
        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)

        if self.ctx.dry_run:
            return inputs

        future = self.submit(PwBaseWorkChain, **inputs)
        self.report(f'launching SCF PwBaseWorkChain<{future.pk}>')

        return ToContext(scf_workchain=future)

    def inspect_scf(self):
        """Verify that the PwBaseWorkChain finished successfully."""

        workchain = self.ctx.scf_workchain
        if not workchain.is_finished_ok:
            self.report(f'SCF PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF

    def should_run_upf2plotcore(self):
        """Don't calculate the core wavefunction data if one has already been provided."""

        return 'core_wfc_data' not in self.inputs

    def should_run_replot(self):
        """Run the WorkChain as a two-step production + replot process if requested."""

        return self.inputs.run_replot.value

    def run_upf2plotcore(self):
        """Generate the core-wavefunction data on-the-fly, if no data is given in the inputs.

        This will determine which pseudopotential is assigned to the atomic species of the
        same element as the absorbing atom, though not the absorbing atom itself, thus the
        corresponding species must use a pseudopotential which contains the correct GIPAW
        information required by the upf2plotcore.sh helper script.

        As this uses the AiiDA-Shell plugin, we assume that this is already installed.
        """

        ShellJob = CalculationFactory('core.shell')  # pylint: disable=invalid-name

        pw_inputs = self.exposed_inputs(PwBaseWorkChain, 'scf')
        pseudo_dict = pw_inputs['pw']['pseudos']
        abs_kind = self.ctx.abs_kind

        upf = pseudo_dict[abs_kind.symbol]

        shell_inputs = {}

        shell_inputs['code'] = self.inputs.upf2plotcore_code
        shell_inputs['nodes'] = {'upf': upf}
        shell_inputs['arguments'] = orm.List(list=['{upf}'])
        shell_inputs['metadata'] = {'call_link_label': 'upf2plotcore'}

        shelljob_node = self.submit(ShellJob, **shell_inputs)
        self.report(f'Launching ShellJob for upf2plotcore.sh<{shelljob_node.pk}>')

        return ToContext(upf2plotcore_node=shelljob_node)

    def inspect_upf2plotcore(self):
        """Check that the output from the upf2plotcore step has yielded a meaningful result.

        This will simply check that the core wavefunction data returned contains at least
        one core state and return an error if this is not the case.
        """

        shelljob_node = self.ctx.upf2plotcore_node
        core_wfc_data = shelljob_node.outputs.stdout
        header_line = shelljob_node.outputs.stdout.get_content()[:40]
        num_core_states = int(header_line.split(' ')[5])
        if num_core_states == 0:
            return self.exit_codes.ERROR_NO_GIPAW_INFO_FOUND
        self.ctx.core_wfc_data = core_wfc_data

    def run_all_xspectra_prod(self):
        """Run an `XspectraBaseWorkChain` for each 3-vector given for epsilon."""

        eps_vectors = self.inputs.eps_vectors.get_list()
        parent_folder = self.ctx.scf_workchain.outputs.remote_folder
        core_wfc_data = self.ctx.core_wfc_data

        calc_number = 0
        for calc_number, vector in enumerate(eps_vectors):
            xspectra_inputs = AttributeDict(self.exposed_inputs(XspectraBaseWorkChain, 'xs_prod'))
            xspectra_parameters = xspectra_inputs.xspectra.parameters.get_dict()

            parent_folder = self.ctx.scf_workchain.outputs.remote_folder
            xspectra_inputs.xspectra.parent_folder = parent_folder
            xspectra_inputs.xspectra.core_wfc_data = core_wfc_data
            xspectra_inputs.metadata.call_link_label = f'xas_{calc_number}_prod'

            for index in [0, 1, 2]:
                xspectra_parameters['INPUT_XSPECTRA'][f'xepsilon({index + 1})'] = vector[index]
            xspectra_inputs.xspectra.parameters = orm.Dict(xspectra_parameters)

            if self.ctx.dry_run:
                return xspectra_inputs

            future_xspectra = self.submit(XspectraBaseWorkChain, **xspectra_inputs)
            self.to_context(xspectra_prod_calculations=append_(future_xspectra))
            self.report(
                f'launching XspectraWorkChain<{future_xspectra.pk}> for epsilon vector {vector}'
                ' (Lanczos production)'
            )

    def inspect_all_xspectra_prod(self):
        """Verify that the `XspectraBaseWorkChain` Lanczos production sub-processes finished successfully."""

        calculations = self.ctx.xspectra_prod_calculations
        unrecoverable_failures = False  # pylint: disable=unused-variable

        for calculation in calculations:
            vector = calculation.outputs.output_parameters.get_dict()['xepsilon']
            if not calculation.is_finished_ok:
                self.report(f'XspectraBaseWorkChain <{vector}>'
                            ' failed with exit status {calculation.exit_status}.')
                unrecoverable_failures = True
            else:
                self.report(f'XspectraBaseWorkChain <{vector}> finished successfully.')
                self.ctx['finished_lanczos'].append(calculation)
        if unrecoverable_failures:
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_XSPECTRA

        self.ctx.all_lanczos_computed = True

    def run_all_xspectra_plot(self):
        """Run an `XspectraBaseWorkChain` for each 3-vector given for epsilon to plot the final spectra.

        This part simply convolutes and plots the spectra from the already-computed Lanczos
        of ``run_all_xspectra_plot``. Only run if requested via ``run_replot`` in the inputs.
        """

        finished_calculations = self.ctx.finished_lanczos

        core_wfc_data = self.ctx.core_wfc_data

        for calc_number, parent_xspectra in enumerate(finished_calculations):
            xspectra_inputs = AttributeDict(self.exposed_inputs(XspectraBaseWorkChain, 'xs_plot'))
            # The epsilon vectors are not needed in the case of a replot, however they
            # will be needed by the Parser at the end
            xspectra_parameters = xspectra_inputs.xspectra.parameters.get_dict()

            parent_output_dict = parent_xspectra.outputs.output_parameters.get_dict()
            parent_calc_job = parent_xspectra.outputs.output_parameters.creator
            eps_vector = parent_output_dict['xepsilon']
            xspectra_parameters['INPUT_XSPECTRA']['xepsilon(1)'] = eps_vector[0]
            xspectra_parameters['INPUT_XSPECTRA']['xepsilon(2)'] = eps_vector[1]
            xspectra_parameters['INPUT_XSPECTRA']['xepsilon(3)'] = eps_vector[2]
            xspectra_inputs.xspectra.parent_folder = parent_xspectra.outputs.remote_folder
            xspectra_inputs.kpoints = parent_calc_job.inputs.kpoints
            xspectra_inputs.xspectra.core_wfc_data = core_wfc_data
            xspectra_inputs.metadata.call_link_label = f'xas_{calc_number}_plot'

            xspectra_inputs.xspectra.parameters = orm.Dict(xspectra_parameters)

            if self.ctx.dry_run:
                return xspectra_inputs

            future_xspectra = self.submit(XspectraBaseWorkChain, **xspectra_inputs)
            self.report(
                f'launching XspectraBaseWorkChain<{future_xspectra.pk}> for epsilon vector {eps_vector} (Replot)'
            )
            self.to_context(xspectra_plot_calculations=append_(future_xspectra))

    def inspect_all_xspectra_plot(self):
        """Verify that the `XspectraBaseWorkChain` re-plot sub-processes finished successfully."""

        calculations = self.ctx.xspectra_plot_calculations

        finished_replots = []
        unrecoverable_failures = False
        for calculation in calculations:
            if not calculation.is_finished_ok:
                self.report(f'XspectraBaseWorkChain failed with exit status {calculation.exit_status}')
                unrecoverable_failures = True
            else:
                finished_replots.append(calculation)
        if unrecoverable_failures:
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_XSPECTRA
        self.ctx.finished_replots = finished_replots

    def results(self):
        """Attach the important output nodes to the outputs of the WorkChain.

        This will collect the SCF and XSpectra output parameters, as well as the
        powder spectrum (if requested)
        """

        xspectra_prod_calcs = self.ctx.finished_lanczos
        if self.inputs.run_replot.value:
            final_calcs = self.ctx.finished_replots
        else:
            final_calcs = self.ctx.finished_lanczos

        eps_powder_vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        basis_vectors_present = False
        for calc in final_calcs:
            out_params = calc.outputs.output_parameters
            in_params = calc.inputs.xspectra.parameters.get_dict()
            eps_vectors = out_params['xepsilon']
            xcoordcrys = out_params['xcoordcrys']
            calc_type = in_params['INPUT_XSPECTRA']['calculation']
            if xcoordcrys is False:
                self.report(
                    'WARNING: calculations were set to use a cartesian basis instead of a '
                    'crystallographic one. Please use ``"xcoordcrys" : True`` to compute the '
                    'powder spectrum using this WorkChain.'
                )
                break
            if eps_vectors in eps_powder_vectors and calc_type == 'xanes_dipole':
                basis_vectors_present = True

        eps_basis_calcs = {}
        if self.inputs.get_powder_spectrum and basis_vectors_present:
            a_vector_present = False
            b_vector_present = False
            c_vector_present = False
            for plot_calc in final_calcs:
                out_params = plot_calc.outputs.output_parameters
                plot_vector = out_params['xepsilon']
                spectrum_node = plot_calc.outputs.spectra
                if plot_vector in eps_powder_vectors:
                    if plot_vector == [1., 0., 0.]:
                        eps_basis_calcs['eps_100'] = spectrum_node
                        a_vector_present = True
                    if plot_vector == [0., 1., 0.]:
                        eps_basis_calcs['eps_010'] = spectrum_node
                        b_vector_present = True
                    if plot_vector == [0., 0., 1.]:
                        eps_basis_calcs['eps_001'] = spectrum_node
                        c_vector_present = True

            # Here, we control for the case where the A and B vectors are given, but C is
            # missing, which would cause a problem for ``get_powder_spectrum``
            if a_vector_present and b_vector_present and not c_vector_present:
                self.report(
                    'WARNING: epsilon vectors for [1.0 0.0 0.0] and [0.0 1.0 0.0] were '
                    'found, but not for [0.0 0.0 1.0]. Please ensure that the vectors '
                    'perpendicular and parallel to the C-axis are defined in the case '
                    'of a system with dichorism.'
                )
            else:
                eps_basis_calcs['metadata'] = {'call_link_label': 'get_powder_spectrum'}
                powder_spectrum = get_powder_spectrum(**eps_basis_calcs)
                self.out('powder_spectrum', powder_spectrum)
        elif self.inputs.get_powder_spectrum and not basis_vectors_present:
            self.report(
                'WARNING: A powder spectrum was requested, but none of the epsilon vectors '
                'given are suitable to compute this.'
            )

        self.out('parameters_scf', self.ctx.scf_workchain.outputs.output_parameters)

        all_xspectra_prod_calcs = {}
        for index, calc in enumerate(xspectra_prod_calcs):
            all_xspectra_prod_calcs[f'xas_{index}'] = calc

        xspectra_prod_params = {}
        for key, node in all_xspectra_prod_calcs.items():
            output_params = node.outputs.output_parameters
            xspectra_prod_params[key] = output_params
        self.out('parameters_xspectra', xspectra_prod_params)

        all_final_spectra = {}
        for index, calc in enumerate(final_calcs):
            all_final_spectra[f'xas_{index}'] = calc.outputs.spectra

        all_final_spectra['metadata'] = {'call_link_label': 'merge_spectra'}
        output_spectra = merge_spectra(**all_final_spectra)

        self.out('spectra', output_spectra)

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
