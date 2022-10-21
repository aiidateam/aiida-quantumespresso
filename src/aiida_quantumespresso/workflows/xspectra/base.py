# -*- coding: utf-8 -*-
"""Workchain to compute the X-ray absorption spectrum for a given structure.

Uses QuantumESPRESSO pw.x and xspectra.x, requires ``AiiDA-Shell`` to run ``upf2plotcore.sh``.
"""
from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, calcfunction, if_, while_
from aiida.orm.nodes.data.base import to_aiida_type
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory

from aiida_quantumespresso.calculations.functions.create_kpoints_from_distance import create_kpoints_from_distance
from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.workflows.protocols.utils import ProtocolMixin, recursive_merge

PwCalculation = CalculationFactory('quantumespresso.pw')
PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
XspectraCalculation = CalculationFactory('quantumespresso.xspectra')
XyData = DataFactory('array.xy')


@calcfunction
def get_all_spectra(**kwargs):
    """Compile all calculated spectra into a single XyData node for easier plotting.

    This will take all the data for each spectrum produced during the re-plot step and output
    them into a single XyData node.
    """

    output_spectra = XyData()
    y_arrays_list = []
    y_units_list = []
    y_labels_list = []

    spectra = [node for label, node in kwargs.items() if label != 'metadata']

    for spectrum_node in spectra:
        calc_node = spectrum_node.creator
        calc_out_params = calc_node.res
        eps_vector = calc_out_params['epsilon_vector']

        old_y_component = spectrum_node.get_y()
        if len(old_y_component) == 1:
            y_array = old_y_component[0][1]
            y_units = old_y_component[0][2]
            y_arrays_list.append(y_array)
            y_units_list.append(y_units)
            y_labels_list.append(f'sigma_{eps_vector[0]}_{eps_vector[1]}_{eps_vector[2]}')
        elif len(old_y_component) == 3:
            y_tot = old_y_component[0][1]
            y_tot_units = old_y_component[0][2]
            y_tot_label = f'sigma_tot_{eps_vector[0]}_{eps_vector[1]}_{eps_vector[2]}'
            y_arrays_list.append(y_tot)
            y_units_list.append(y_tot_units)
            y_labels_list.append(y_tot_label)

            y_up = old_y_component[1][1]
            y_up_units = old_y_component[1][2]
            y_up_label = f'sigma_up_{eps_vector[0]}_{eps_vector[1]}_{eps_vector[2]}'
            y_arrays_list.append(y_up)
            y_units_list.append(y_up_units)
            y_labels_list.append(y_up_label)

            y_down = old_y_component[2][1]
            y_down_units = old_y_component[2][2]
            y_down_label = f'sigma_down_{eps_vector[0]}_{eps_vector[1]}_{eps_vector[2]}'
            y_arrays_list.append(y_down)
            y_units_list.append(y_down_units)
            y_labels_list.append(y_down_label)

        x_array = spectrum_node.get_x()[1]
        x_label = spectrum_node.get_x()[0]
        x_units = spectrum_node.get_x()[2]

    output_spectra.set_x(x_array, x_label, x_units)
    output_spectra.set_y(y_arrays_list, y_labels_list, y_units_list)

    return output_spectra


@calcfunction
def get_powder_spectrum(**kwargs):
    """Combine the output spectra into a single "Powder" spectrum, representing the K-edge XAS of a powder sample.

    Note that this step should only be requested (``inputs.get_powder_spectrum= True``, defaults to
    ``False``) for valid crystal structures. The function interprets the crystal symmetry to be
    exploited based on the polarisation (epsilon) vectors given to the function, since (by
    default) the WorkChain calls for calculations of all three basis vectors and thus it is
    assumed that the user would know the symmetry of the structure before deciding to set the
    vectors manually. If the polarisation vectors supplied to the ``WorkChain`` would cause a
    failure in the ``CalcFunction`` step, then the step is skipped entirely with a warning (see
    ``cls.results``)
    """

    spectra = {label: node for label, node in kwargs.items() if label != 'metadata'}

    # If the system is isochoric (e.g. a cubic system) then the three basis vectors are equal
    # to each other, thus we simply return the
    if len(spectra) == 1:
        vectors = list(spectra.keys())
        powder_spectrum = spectra[vectors[0]]
        powder_x = powder_spectrum.get_x()[1]
        powder_y = powder_spectrum.get_y()[0][1]

        powder_data = orm.XyData()
        powder_data.set_x(powder_x, 'energy', 'eV')
        powder_data.set_y(powder_y, 'sigma', 'n/a')

    # if the system is dichoric (e.g. a hexagonal system) then the A and B periodic
    # dimensions are equal to each other by symmetry, thus the powder spectrum is simply
    # the average of 2x the 1 0 0 eps vector and 1x the 0 0 1 eps vector
    if len(spectra) == 2:
        # Since the individual vectors are labelled, we can extract just the spectra needed
        # to produce the powder and leave the rest
        vectors = list(spectra.keys())

        for label in vectors:
            if label in ['eps_100', 'eps_010']:
                spectrum_a = spectra[label]
            elif label in ['eps_001']:
                spectrum_c = spectra[label]

        powder_x = spectrum_a.get_x()[1]
        yvals_a = spectrum_a.get_y()[0][1]
        yvals_c = spectrum_c.get_y()[0][1]

        powder_y = ((yvals_a * 2) + yvals_c) / 3
        powder_data = orm.XyData()
        powder_data.set_x(powder_x, 'energy', 'eV')
        powder_data.set_y(powder_y, 'sigma', 'n/a')

    # if the system is trichoric (e.g. a monoclinic system) then no periodic dimensions
    # are equal by symmetry, thus the powder spectrum is the average of the three basis
    # dipole vectors (1 0 0, 0 1 0, 0 0 1)
    if len(spectra) == 3:
        # Since the individual vectors are labelled, we can extract just the spectra needed to
        # produce the powder and leave the rest
        spectrum_a = spectra['eps_100']
        spectrum_b = spectra['eps_010']
        spectrum_c = spectra['eps_001']

        powder_x = spectrum_a.get_x()[1]
        yvals_a = spectrum_a.get_y()[0][1]
        yvals_b = spectrum_b.get_y()[0][1]
        yvals_c = spectrum_c.get_y()[0][1]

        powder_y = (yvals_a + yvals_b + yvals_c) / 3

        powder_data = orm.XyData()
        powder_data.set_x(powder_x, 'energy', 'eV')
        powder_data.set_y(powder_y, 'sigma', 'n/a')

    return powder_data


def validate_scf(inputs, _):
    """Validate the scf parameters."""
    parameters = inputs['pw']['parameters'].get_dict()
    if parameters.get('CONTROL', {}).get('calculation', 'scf') != 'scf':
        return '`CONTOL.calculation` in `scf.pw.parameters` is not set to `scf`.'


def validate_inputs(inputs, _):
    """Validate the inputs before launching the WorkChain."""

    eps_vector_list = inputs['eps_vectors'].get_list()
    if len(eps_vector_list) == 0:
        return 'Error: eps_vectors list empty.'
    if 'core_wfc_data' not in inputs and 'upf2plotcore_code' not in inputs:
        return 'Error: either a core wavefunction file or a code node for upf2plotcore.sh must be provided.'
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


class XspectraBaseWorkChain(ProtocolMixin, WorkChain):
    """Workchain to compute X-ray absorption spectra for a given structure using Quantum ESPRESSO.

    The workflow follows the process required to compute the XAS of an input structure:
    an SCF calculation is performed using the provided structure, which is then followed by
    the calculation of the XAS itself by XSpectra. The calculations performed by the WorkChain
    in a typical run will be:
        * PwSCF calculation with pw.x of the input structure with a core-hole present.
        * Generation of core-wavefunction data with upf2plotcore.sh (if requested).
        * XAS calculation with xspectra.x to compute the Lanczos coefficients at a high level
          of broadening (0.8 eV by default).
        * XAS calculation with xspectra.x to re-plot the XANES spectrum from the computed
          Lanczos coefficients at a lower level of broadening (0.3 eV by default).
        * Collation of output data from pw.x and xspectra.x calculations, including a
          combination of XANES dipole spectra based on polarisation vectors to represent the
          powder spectrum of the structure (if requested).

    The radial part of the core wavefunction (i.e. the atomic state containing the core-hole
    in the absorbing atom) derived from the ground-state of the absorbing element can be
    provided as a top-level input or produced by the WorkChain. If left to the WorkChain,
    the ground-state pseudopotential assigned to the absorbing element will be used to generate
    this data using the upf2plotcore.sh utility script (via the ``AiiDA-Shell`` plugin).

    In its current stage of development, the workflow requires the following:
        * An input structure where the desired absorbing atom in the system is marked in
          some way. The default behaviour for the WorkChain is to set this as 'X', however
          this can be changed via the `overrides` dictionary.
        * A code node for ``upf2plotcore``, configured for the ``AiiDA-Shell`` plugin
          (https://github.com/sphuber/aiida-shell). Alternatively, a ``SinglefileData`` node
          from a previous ``ShellJob`` run can be supplied under ``builder.core_wfc_data``.
        * A suitable pair of pseudopotentials for the element type of the absorbing atom,
          one for the ground-state occupancy which contains GIPAW informtation for the core
          level of interest for the XAS (e.g. 1s in the case of a K-edge calculation) and
          the other containing a core hole.
            * For the moment this can be passed either via the ``core_hole_pseudos`` field
              in ``get_builder_from_protocol`` or via the overrides, but will be
              changed later once full families of core-hole pseudopotentials become
              available.
    """

    # pylint: disable=too-many-public-methods, too-many-statements

    @classmethod
    def define(cls, spec):
        """Define the process specification."""

        super().define(spec)
        # yapf: disable
        spec.expose_inputs(
            PwBaseWorkChain,
            namespace='scf',
            exclude=('pw.parent_folder', 'pw.structure', 'clean_workdir'),
            namespace_options={
                'help': ('Input parameters for the `pw.x` calculation.'),
                'validator': validate_scf,
            }
        )
        spec.expose_inputs(
            XspectraCalculation,
            namespace='xs_prod',
            exclude=('parent_folder', 'kpoints', 'core_wfc_data'),
            namespace_options={'help': ('Input parameters for the initial `xspectra.x` calculation'
           ' to compute the Lanczos.')}
        )
        spec.expose_inputs(
            XspectraCalculation,
            namespace='xs_plot',
            exclude=('parent_folder', 'kpoints', 'core_wfc_data'),
            namespace_options={'help': ('Input parameters for the re-plot `xspectra.x` calculation of the Lanczos.')}
        )
        spec.inputs.validator = validate_inputs
        spec.input(
            'structure',
            valid_type=orm.StructureData,
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
            'upf2plotcore_code',
            valid_type=orm.Code,
            required=False,
            help='The code node required for upf2plotcore.sh configured for ``AiiDA-Shell``. '
            'Must be provided if `core_wfc_data` is not provided.'
        )
        spec.input(
            'clean_workdir',
            valid_type=orm.Bool,
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
        spec.input('kpoints', valid_type=orm.KpointsData, required=False,
            help='An explicit k-points list or mesh for the xspectra.x calculations. Either this or `kpoints_distance`'
                   ' has to be provided.')
        spec.input('kpoints_distance', valid_type=orm.Float, required=False,
            help='The minimum desired distance in 1/Å between k-points in reciprocal space. The explicit k-points will '
                 'be generated automatically by a calculation function based on the input structure. '
                   'Applies only to the xspectra.x calculations.')
        spec.input('kpoints_force_parity', valid_type=orm.Bool, required=False,
            help='Optional input when constructing the k-points based on a desired `kpoints_distance`. Setting this to '
                 '`True` will force the k-point mesh to have an even number of points along each lattice vector except '
                 'for any non-periodic directions. Applies only to the xspectra.x calculations.')
        spec.outline(
            cls.setup,
            cls.validate_kpoints,
            cls.run_scf,
            cls.inspect_scf,
            if_(cls.should_run_upf2plotcore)(cls.run_upf2plotcore),
            while_(cls.should_repeat_xs_prod)(
                cls.run_all_xspectra_prod,
                cls.inspect_all_xspectra_prod,
            ),
            cls.run_all_xspectra_plot,
            cls.inspect_all_xspectra_plot,
            cls.results,
        )

        spec.exit_code(202, 'ERROR_INVALID_INPUT_KPOINTS',
            message='Neither the `kpoints` nor the `kpoints_distance` input was specified.')
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_SCF', message='The SCF sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_XSPECTRA', message='One or more XSpectra sub processes failed')
        spec.output('output_parameters_scf', valid_type=orm.Dict, help='The output parameters of the SCF'
                    ' `PwBaseWorkChain`.')
        spec.output_namespace(
            'output_parameters_xspectra',
            valid_type=orm.Dict,
            help='The output dictionaries of each `XspectraCalculation` performed',
            dynamic=True
        )
        spec.output(
            'output_spectra',
            valid_type=orm.XyData,
            help='An XyData node containing all the final spectra produced by the WorkChain.'
        )
        spec.output('powder_spectrum', valid_type=orm.XyData, required=False, help='The simulated powder spectrum')
        # yapf: disable
    @classmethod
    def get_protocol_filepath(cls):
        """Return ``pathlib.Path`` to the ``.yaml`` file that defines the protocols."""
        from importlib_resources import files

        # pylint: disable=relative-beyond-top-level
        from ..protocols import xspectra as xs_protocols
        return files(xs_protocols) / 'base.yaml'
# pylint: enable=relative-beyond-top-level
    @classmethod
    def get_builder_from_protocol(
        cls, pw_code, xs_code, structure, core_wfc_data=None, upf2plotcore_code=None,
        core_hole_pseudos=None, pw_protocol=None, xs_protocol=None, overrides=None,
        options=None, **kwargs
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param pw_code: the ``Code`` instance configured for the ``quantumespresso.pw``
                        plugin.
        :param xs_code: the ``Code`` instance configured for the
                        ``quantumespresso.xspectra`` plugin.
        :param upf2plotcore_code: the AiiDA-Shell ``Code`` instance configured for the
                                  upf2plotcore.sh shell script, used to generate the core
                                  wavefunction data.
        :param structure: the ``StructureData`` instance to use.
        :param core_hole_pseudos: the core-hole pseudopotential pair (ground-state and
                                  excited-state) for the chosen absorbing atom.
        :param pw_protocol: the px.x protocol to use. If not specified, the default set by
                            the PwBaseWorkChain will be used instead.
        :param overrides: optional dictionary of inputs to override the defaults of the
                          XspectraWorkChain itself.
        :param xs_protocol: the xspectra.x protocol to use, which defines pw.x settings
                            to be used for the core-hole treatment. Defaults to the
                            "Full-core-hole" approach ("full") if not specified.
        :param kwargs: additional keyword arguments that will be passed to the
            ``get_builder_from_protocol`` of all the sub processes that are called by this
            workchain.
        :return: a process builder instance with all inputs defined ready for launch.
        """

        # Get the default inputs from the PwBaseWorkChain and override them with those
        # required for the chosen core-hole treatment
        inputs = recursive_merge(
            left=PwBaseWorkChain.get_protocol_inputs(protocol=pw_protocol),
            right=cls.get_protocol_inputs(protocol=xs_protocol, overrides=overrides)
        )

        args = (pw_code, structure, pw_protocol)
        scf = PwBaseWorkChain.get_builder_from_protocol(
            *args, overrides=inputs.get('scf', None), options=options, **kwargs
        )

        scf.pop('clean_workdir', None)
        scf['pw'].pop('structure', None)
        abs_atom_marker = inputs['abs_atom_marker']
        for kind in structure.kinds:
            if kind.name == abs_atom_marker:
                abs_element = kind.symbol

        # pylint: disable=no-member
        builder = cls.get_builder()
        builder.scf = scf
        metadata_xs_prod = inputs.get('xs_prod', {}).get('metadata', {'options': {}})
        metadata_xs_plot = inputs.get('xs_plot', {}).get('metadata', {'options': {}})

        if options:
            metadata_xs_prod['options'] = recursive_merge(metadata_xs_prod['options'], options)
            metadata_xs_plot['options'] = recursive_merge(metadata_xs_plot['options'], options)
        if core_hole_pseudos:
            builder.scf.pw.pseudos[abs_atom_marker] = core_hole_pseudos[abs_atom_marker]
            builder.scf.pw.pseudos[abs_element] = core_hole_pseudos[abs_element]

        builder.xs_prod.code = xs_code
        builder.xs_prod.parameters = orm.Dict(inputs.get('xs_prod', {}).get('parameters'))
        builder.xs_prod.metadata = metadata_xs_prod
        builder.xs_plot.code = xs_code
        builder.xs_plot.parameters = orm.Dict(inputs.get('xs_plot', {}).get('parameters'))
        builder.xs_plot.metadata = metadata_xs_plot

        if upf2plotcore_code:
            builder.upf2plotcore_code = upf2plotcore_code
        elif core_wfc_data:
            builder.core_wfc_data = core_wfc_data

        builder.structure = structure
        builder.eps_vectors = orm.List(list=inputs['eps_vectors'])
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])
        builder.get_powder_spectrum = orm.Bool(inputs['get_powder'])
        builder.abs_atom_marker = orm.Str(inputs['abs_atom_marker'])
        builder.kpoints_distance = orm.Float(inputs['kpoints_distance'])
        builder.kpoints_force_parity = orm.Bool(inputs['kpoints_force_parity'])
        # pylint: enable=no-member
        return builder

    def setup(self):
        """Initialize context variables that are used during the logical flow of the workchain."""

        # self.ctx.serial_clean = 'serial_clean' in self.inputs and self.inputs.serial_clean.value
        self.ctx.dry_run = 'dry_run' in self.inputs and self.inputs.dry_run.value
        self.ctx.all_lanczos_computed = False
        self.ctx.lanczos_to_restart = []
        self.ctx.finished_lanczos = []

    def validate_kpoints(self):
        """Validate the inputs related to k-points for the ``XspectraCalculation``s.

        Either an explicit ``KpointsData`` with given mesh/path, or a desired k-points distance should be specified. In
        the case of the latter, the ``KpointsData`` will be constructed for the input ``StructureData`` using the
        ``create_kpoints_from_distance`` calculation function.
        """
        if all(key not in self.inputs for key in ['kpoints', 'kpoints_distance']):
            return self.exit_codes.ERROR_INVALID_INPUT_KPOINTS

        try:
            kpoints = self.inputs.kpoints
        except AttributeError:
            inputs = {
                'structure': self.inputs.structure,
                'distance': self.inputs.kpoints_distance,
                'force_parity': self.inputs.get('kpoints_force_parity', orm.Bool(False)),
                'metadata': {
                    'call_link_label': 'create_kpoints_from_distance'
                }
            }
            kpoints = create_kpoints_from_distance(**inputs)  # pylint: disable=unexpected-keyword-arg

        self.ctx.xspectra_kpoints = kpoints

    def should_repeat_xs_prod(self):
        """Return whether the Lanczos production step is finished or not."""

        return not self.ctx.all_lanczos_computed

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

        # self.ctx.scf_parent_folder = workchain.outputs.remote_folder

    def should_run_upf2plotcore(self):
        """Don't calculate the core wavefunction data if one has already been provided."""

        return 'core_wfc_data' not in self.inputs

    def run_upf2plotcore(self):
        """If no core wavefunction data node is provided, then generate one on-the-fly.

        This will determine which pseudopotential is assigned to the atomic species of the
        same element as the absorbing atom, though not the absorbing atom itself, thus the
        corresponding species must use a pseudopotential which contains the correct GIPAW
        information required by the upf2plotcore.sh helper script.

        As this uses the AiiDA-Shell plugin, we assume that this is already installed.
        """

        ShellJob = CalculationFactory('core.shell') # pylint: disable=invalid-name

        pw_inputs = self.exposed_inputs(PwBaseWorkChain, 'scf')
        absorbing_species = self.inputs.abs_atom_marker.value

        pseudo_dict = pw_inputs['pw']['pseudos']
        upf = pseudo_dict[absorbing_species]

        shell_inputs = {}

        shell_inputs['code'] = self.inputs.upf2plotcore_code
        shell_inputs['files'] = {'upf': upf}
        shell_inputs['arguments'] = orm.List(list=['{upf}'])
        shell_inputs['metadata'] = {'call_link_label': 'upf2plotcore'}

        shelljob_node = self.submit(ShellJob, **shell_inputs)
        self.report(f'Launching ShellJob for upf2plotcore.sh<{shelljob_node.pk}>')

        return ToContext(upf2plotcore_node=shelljob_node)

    def run_all_xspectra_prod(self):
        """Run an `XspectraCalculation` for each 3-vector given for epsilon to produce the Lanczos coefficients."""

        eps_vectors = self.inputs.eps_vectors.get_list()
        structure = self.inputs.structure
        kpoints = self.ctx.xspectra_kpoints
        if 'core_wfc_data' in self.inputs:
            core_wfc_data = self.inputs.core_wfc_data
        else:
            core_wfc_data = self.ctx.upf2plotcore_node.outputs.stdout

        kinds_list = [kind.name for kind in structure.kinds]
        kinds_list.sort()

        # This will work so long as the PwCalculation part adds the atomic species to
        # the input file in alphabetical order
        kind_counter = 1
        for kind in kinds_list:
            if kind == self.inputs.abs_atom_marker.value:
                xiabs = kind_counter
            else:
                kind_counter += 1

        xspectra_prod_calcs = {}
        if 'xspectra_calc_labels' not in self.ctx.keys():
            xspectra_calc_labels = []
        calc_number = 0
        if len(self.ctx.lanczos_to_restart) == 0: # No restarts have been ordered yet
            for vector in eps_vectors:
                xspectra_inputs = AttributeDict(self.exposed_inputs(XspectraCalculation, 'xs_prod'))
                xspectra_parameters = xspectra_inputs.parameters.get_dict()
                xspectra_parameters['INPUT_XSPECTRA']['xiabs'] = xiabs
                max_walltime_seconds = self.inputs.xs_prod.metadata.options.max_wallclock_seconds
                xspectra_parameters['INPUT_XSPECTRA']['time_limit'] = max_walltime_seconds * 0.9

                xspectra_inputs.parent_folder = self.ctx.scf_workchain.outputs.remote_folder
                xspectra_inputs.kpoints = kpoints
                xspectra_inputs.core_wfc_data = core_wfc_data
                label = f'xas_{calc_number}'
                xspectra_inputs.metadata.call_link_label = f'{label}_prod_iter_1'
                xspectra_inputs.metadata.label = f'{label}_prod_iter_1'
                xspectra_calc_labels.append(label)

                for index in [0, 1, 2]:
                    xspectra_parameters['INPUT_XSPECTRA'][f'xepsilon({index + 1})'] = vector[index]
                xspectra_inputs.parameters = orm.Dict(dict=xspectra_parameters)

                if self.ctx.dry_run:
                    return xspectra_inputs

                future_xspectra = self.submit(XspectraCalculation, **xspectra_inputs)
                self.report(
                    f'launching XspectraCalculation<{future_xspectra.pk}> for epsilon vector {vector}'
                    ' (Lanczos production) (iteration #1)'
                )
                xspectra_prod_calcs[f'{label}_prod'] = future_xspectra
                calc_number += 1

        else: # Some calculations need restarting, so we process these instead
            for calculation in self.ctx.lanczos_to_restart:
                xspectra_inputs = AttributeDict(self.exposed_inputs(XspectraCalculation, 'xs_prod'))
                xspectra_parameters = calculation.inputs.parameters.get_dict()
                xspectra_parameters['INPUT_XSPECTRA']['restart_mode'] = 'restart'
                vector = [
                    xspectra_parameters['INPUT_XSPECTRA']['xepsilon(1)'],
                    xspectra_parameters['INPUT_XSPECTRA']['xepsilon(2)'],
                    xspectra_parameters['INPUT_XSPECTRA']['xepsilon(3)']
                ]

                xspectra_inputs.parent_folder = calculation.outputs.remote_folder
                xspectra_inputs.kpoints = kpoints
                xspectra_inputs.core_wfc_data = core_wfc_data
                parent_label_pieces = calculation.label.split('_')
                iteration = int(parent_label_pieces[-1])
                iteration += 1
                label = f'{parent_label_pieces[0]}_{parent_label_pieces[1]}_{parent_label_pieces[2]}'
                xspectra_inputs.metadata.call_link_label = f'{label}_iter_{iteration}'
                xspectra_inputs.metadata.label = f'{label}_iter_{iteration}'

                xspectra_inputs.parameters = orm.Dict(dict=xspectra_parameters)

                if self.ctx.dry_run:
                    return xspectra_inputs

                future_xspectra = self.submit(XspectraCalculation, **xspectra_inputs)
                self.report(
                    f'launching XspectraCalculation<{future_xspectra.pk}> for epsilon vector {vector}'
                    f' (Lanczos production) (iteration #{iteration})'
                )
                xspectra_prod_calcs[f'{label}'] = future_xspectra

        if 'xspectra_calc_labels' not in self.ctx.keys():
            self.ctx.xspectra_calc_labels = xspectra_calc_labels
        return ToContext(**xspectra_prod_calcs)

    def inspect_all_xspectra_prod(self):
        """Verify that the `XspectraCalculation` Lanczos production sub-processes finished successfully."""

        calculations = []
        labels = self.ctx.xspectra_calc_labels
        for label in labels:
            calculation = self.ctx[f'{label}_prod']
            calculations.append(calculation)
        unrecoverable_failures = False # pylint: disable=unused-variable
        restarts_required = False

        lanczos_to_restart = []
        for calculation in calculations:
            if not calculation.is_finished_ok:
                # Check if any calculations failed due to out-of-walltime
                label_pieces = calculation.label.split('_')
                report_label = f'{label_pieces[0]}_{label_pieces[1]}_{label_pieces[2]}'
                if calculation.exit_status == 400:
                    lanczos_to_restart.append(calculation)
                    self.report(f'Lanczos calculation {report_label} added to restart list.')
                    restarts_required = True
                else:
                    self.report(f'XspectraCalculation <{report_label}>'
                                ' failed with exit status {calculation.exit_status}.')
                    unrecoverable_failures = True
                    return self.exit_codes.ERROR_SUB_PROCESS_FAILED_XSPECTRA
            else:
                self.ctx['finished_lanczos'].append(calculation)
        if restarts_required: # If there are restarts, then pass them to the context
            self.ctx['lanczos_to_restart'] = lanczos_to_restart
        else: # If all is well, then we can continue to the replot stage
            self.ctx.all_lanczos_computed = True

    def run_all_xspectra_plot(self):
        """Run an `XspectraCalculation` for each 3-vector given for epsilon to plot the final spectra.

        This part simply prints the spectra from the already-computed Lanczos stored in.
        """

        finished_calculations = self.ctx.finished_lanczos
        kpoints = self.ctx.xspectra_kpoints

        if 'core_wfc_data' in self.inputs:
            core_wfc_data = self.inputs.core_wfc_data
        else:
            core_wfc_data = self.ctx.upf2plotcore_node.outputs.stdout

        xspectra_plot_calcs = {}
        for parent_xspectra in finished_calculations:
            xspectra_inputs = AttributeDict(self.exposed_inputs(XspectraCalculation, 'xs_plot'))
            # The epsilon or k vectors are not needed in the case of a replot, however they
            # will be needed by the Parser at the end
            xspectra_parameters = xspectra_inputs.parameters.get_dict()

            parent_output_dict = parent_xspectra.res
            eps_vector = parent_output_dict['epsilon_vector']
            xspectra_parameters['INPUT_XSPECTRA']['xepsilon(1)'] = eps_vector[0]
            xspectra_parameters['INPUT_XSPECTRA']['xepsilon(2)'] = eps_vector[1]
            xspectra_parameters['INPUT_XSPECTRA']['xepsilon(3)'] = eps_vector[2]
            parent_label_pieces = parent_xspectra.label.split('_')
            old_label = f'{parent_label_pieces[0]}_{parent_label_pieces[1]}_{parent_label_pieces[2]}'
            new_label = old_label.replace('prod', 'plot')
            xspectra_inputs.parent_folder = parent_xspectra.outputs.remote_folder
            xspectra_inputs.kpoints = kpoints
            xspectra_inputs.core_wfc_data = core_wfc_data
            xspectra_inputs.metadata.label = new_label
            xspectra_inputs.metadata.call_link_label = new_label

            xspectra_inputs.parameters = orm.Dict(dict=xspectra_parameters)

            if self.ctx.dry_run:
                return xspectra_inputs

            future_xspectra = self.submit(XspectraCalculation, **xspectra_inputs)
            self.report(f'launching XspectraCalculation<{future_xspectra.pk}> for epsilon vector {eps_vector} (Replot)')
            xspectra_plot_calcs[f'{new_label}'] = future_xspectra

        return ToContext(**xspectra_plot_calcs)

    def inspect_all_xspectra_plot(self):
        """Verify that the `XspectraCalculation` re-plot sub-processes finished successfully."""

        calculations = []
        labels = self.ctx.xspectra_calc_labels
        for label in labels:
            calculation = self.ctx[f'{label}_plot']
            calculations.append(calculation)
        unrecoverable_failures = False

        finished_replots = []
        for calculation in calculations:
            if not calculation.is_finished_ok:
                self.report(f'XspectraCalculation failed with exit status {calculation.exit_status}')
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
        xspectra_plot_calcs = self.ctx.finished_replots

        eps_powder_vectors = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
        for calc in xspectra_plot_calcs:
            out_params = calc.res
            in_params = calc.inputs.parameters.get_dict()
            eps_vectors = out_params['epsilon_vector']
            coord_system = out_params['vector_coord_system']
            calc_type = in_params['INPUT_XSPECTRA']['calculation']
            if coord_system == 'cartesian':
                self.report(
                    'WARNING: calculations were set to use a cartesian basis instead of a '
                    'crystallographic one. Please use ``"xcoordcrys" : True`` to compute the '
                    'powder spectrum.'
                )
                basis_vectors_present = False
                break
            if eps_vectors in eps_powder_vectors and calc_type == 'xanes_dipole':
                basis_vectors_present = True

        eps_basis_calcs = {}
        if self.inputs.get_powder_spectrum and basis_vectors_present:
            a_vector_present = False
            b_vector_present = False
            c_vector_present = False
            for plot_calc in xspectra_plot_calcs:
                out_params = plot_calc.res
                plot_vector = out_params['epsilon_vector']
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
                eps_basis_calcs['metadata'] = {'call_link_label': 'compile_powder'}
                powder_spectrum = get_powder_spectrum(**eps_basis_calcs)
                self.out('powder_spectrum', powder_spectrum)
        elif self.inputs.get_powder_spectrum and not basis_vectors_present:
            self.report(
                'WARNING: A powder spectrum was requested, but none of the epsilon vectors '
                'given are suitable to compute this.'
            )

        self.out('output_parameters_scf', self.ctx.scf_workchain.outputs.output_parameters)

        all_xspectra_prod_calcs = {}
        for calc in xspectra_prod_calcs:
            label_pieces = calc.label.split('_')[:-2]
            label = f'{label_pieces[0]}_{label_pieces[1]}_{label_pieces[2]}'
            all_xspectra_prod_calcs[label] = calc

        xspectra_prod_params = {}
        for key, node  in all_xspectra_prod_calcs.items():
            output_params = node.outputs.output_parameters
            xspectra_prod_params[key] = output_params
        self.out('output_parameters_xspectra', xspectra_prod_params)

        all_final_spectra = {}
        for calc in xspectra_plot_calcs:
            all_final_spectra[calc.label] = calc.outputs.spectra

        all_final_spectra['metadata'] = {'call_link_label': 'compile_all_spectra'}
        output_spectra = get_all_spectra(**all_final_spectra)

        self.out('output_spectra', output_spectra)
        if self.inputs.clean_workdir.value is True:
            self.report('workchain succesfully completed, cleaning remote folders')
        elif self.inputs.clean_workdir.value is False:
            self.report('workchain succesfully completed, remote folders will be kept')

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
