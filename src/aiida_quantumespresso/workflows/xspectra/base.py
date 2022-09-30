# -*- coding: utf-8 -*-
"""Workchain to compute the X-ray absorption spectrum for a given structure.

Uses QuantumESPRESSO pw.x and xspectra.x.
"""
from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import ToContext, WorkChain, calcfunction, if_, while_
from aiida.plugins import CalculationFactory, DataFactory, WorkflowFactory

from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.workflows.protocols.utils import ProtocolMixin, recursive_merge

PwCalculation = CalculationFactory('quantumespresso.pw')
PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
XspectraCalculation = CalculationFactory('quantumespresso.xspectra')
XyData = DataFactory('array.xy')


@calcfunction
def get_all_spectra(**kwargs):
    """Compile all calculated spectra into a single XyData node for easier plotting.

    This will take only the "total sigma" value for each spectrum produced during the
    re-plot step and output a single XyData node containing all the obtained spectra.
    """

    output_spectra = XyData()
    y_arrays_list = []
    y_units_list = []
    y_labels_list = []

    calculations = [orm.load_node(pk) for label, pk in kwargs.items() if label != 'metadata']

    for calc in calculations:
        spectrum_node = calc.get_outgoing().get_node_by_label('spectra')
        parent_calc = calc.get_incoming().get_node_by_label('parent_folder').creator
        parent_out_params = parent_calc.get_outgoing().get_node_by_label('output_parameters').get_dict()
        eps_vector = parent_out_params['epsilon_vector']

        old_y_component = spectrum_node.get_y()[0]
        y_array = old_y_component[1]
        y_units = old_y_component[2]
        y_arrays_list.append(y_array)
        y_units_list.append(y_units)
        y_labels_list.append(f'sigma_{eps_vector[0]}_{eps_vector[1]}_{eps_vector[2]}')

    x_array = spectrum_node.get_x()[1]
    x_label = spectrum_node.get_x()[0]
    x_units = spectrum_node.get_x()[2]
    output_spectra.set_x(x_array, x_label, x_units)
    output_spectra.set_y(y_arrays_list, y_labels_list, y_units_list)

    return output_spectra


@calcfunction
def get_powder_spectrum(**kwargs):
    """Combine the output spectra into a single "Powder" spectrum.

    Note that this step should only be requested (``get_powder=True``) if the chosen
    crystal system meaningfully obeys the rules for either isochorism, dichorism, or
    trichorism.
    """

    calculations = {label: orm.load_node(pk) for label, pk in kwargs.items() if label != 'metadata'}

    # If the system is isochoric (e.g. a cubic system) then the three basis vectors are equal
    # to each other, thus we simply return the
    if len(calculations) == 1:
        vectors = list(calculations.keys())
        powder_spectrum = calculations[vectors[0]].get_outgoing().get_node_by_label('spectra')
        powder_x = powder_spectrum.get_x()[1]
        powder_y = powder_spectrum.get_y()[0][1]

        powder_data = orm.XyData()
        powder_data.set_x(powder_x, 'energy', 'eV')
        powder_data.set_y(powder_y, 'sigma', 'n/a')

    # if the system is dichoric (e.g. a hexagonal system) then the A and B periodic
    # dimensions are equal to each other by symmetry, thus the powder spectrum is simply
    # the average of 2x the 1 0 0 eps vector and 1x the 0 0 1 eps vector
    if len(calculations) == 2:
        # Since the individual vectors are labelled, we can extract just the spectra needed to
        # produce the powder and leave the rest
        vectors = list(calculations.keys())

        for label in vectors:
            if label in ['eps_100', 'eps_010']:
                spectrum_a = calculations[label].get_outgoing().get_node_by_label('spectra')
            elif label in ['eps_001']:
                spectrum_c = calculations[label].get_outgoing().get_node_by_label('spectra')

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
    if len(calculations) == 3:
        # Since the individual vectors are labelled, we can extract just the spectra needed to
        # produce the powder and leave the rest
        spectrum_a = calculations['eps_100'].get_outgoing().get_node_by_label('spectra')
        spectrum_b = calculations['eps_010'].get_outgoing().get_node_by_label('spectra')
        spectrum_c = calculations['eps_001'].get_outgoing().get_node_by_label('spectra')

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
        return 'Error: either a core wavefunction file or a code node for upf2plotcore must be provided.'


class XspectraBaseWorkChain(ProtocolMixin, WorkChain):
    """Workchain to compute X-ray absorption spectra for a given structure using Quantum ESPRESSO.

    The workflow follows the process required to compute the XAS of an input structure:
    an SCF calculation is performed using the provided structure, which is then followed by
    the calculation of the XAS itself by XSpectra.

    The radial part of the core wavefunction (i.e. the atomic state containing the core-hole
    in the absorbing atom) derived from the ground-state of the absorbing element can be
    provided as a top-level input or produced by the WorkChain. If left to the WorkChain, the
    ground-state pseudopotential assigned to the absorbing element will be used to generate
    this data using the upf2plotcore.sh utility script (via the AiiDA-Shell plugin).

    In its current stage of development, the workflow requires the following:
        * An input structure where the desired absorbing atom in the system is marked in some
          way. The default behaviour for the WorkChain is to set this as 'X' by default,
          however this can be changed via the `overrides` dictionary.
        * An already-installed code node for 'upf2plotcore.sh@localhost', produced by running
          the AiiDA-Shell plugin with upf2plotcore.sh once via the launch_shell_job() method.
          Alternatively, a pre-computed `stdout` node can be supplied under
          `builder.core_wfc_data`.
        * A suitable pair of pseudopotentials for the element type of the absorbing atom,
          one for the ground-state occupancy which contains GIPAW informtation on the core
          level of interest for the XAS (e.g. 1s in the case of a K-edge calculation) and the
          other containing a core hole.
            * For the moment, this must be passed through the use of overrides in the
              get_builder_from_protocol step, but will be changed later once full families of
              core-hole pseudopotentials become available.
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
            namespace_options={'help': ('Input parameters for the `xspectra.x` calculation.')}
        )
        spec.expose_inputs(
            XspectraCalculation,
            namespace='xs_plot',
            exclude=('parent_folder', 'kpoints', 'core_wfc_data'),
            namespace_options={'help': ('Input parameters for the `xspectra.x` calculation.')}
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
            default=lambda: orm.List(list=[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
            help=(
                'The list of 3-vectors to use in XSpectra sub-processes. '
                'The number of sub-lists will subsequently define the number of XSpectra calculations to perform'
            ),
        )
        spec.input(
            'clean_workdir',
            valid_type=orm.Bool,
            default=lambda: orm.Bool(False),
            help=('If `True`, work directories of all called calculation will be cleaned at the end of execution.'),
        )
        spec.input(
            'abs_atom_marker',
            valid_type=orm.Str,
            default=lambda: orm.Str('X'),
            help=(
                'The name for the Kind representing the absorbing atom in the structure. '
                'Must corespond to a Kind within the StructureData node supplied to the calculation.'
            ),
        )
        spec.input(
            'collect_powder',
            valid_type=orm.Bool,
            default=lambda: orm.Bool(False),
            help=(
                'If `True`, the WorkChain will combine XAS spectra computed using the XAS basis vectors defined '
                'according to the `get_powder_spectrum` CalcFunction.'
            ),
        )
        spec.input(
            'core_wfc_data',
            valid_type=orm.SinglefileData,
            required=False,
            help='The core wavefunction data file extracted from the ground-state pseudo for the absorbing atom'
        )
        spec.input(
            'upf2plotcore_code',
            valid_type=orm.Code,
            required=False,
            help='The code node required for upf2plotcore.sh. Must be provided if `core_wfc_data` is not provided'
        )
        spec.outline(
            cls.setup,
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

        from ..protocols import xspectra as xs_protocols  # pylint: disable=relative-beyond-top-level
        return files(xs_protocols) / 'base.yaml'

    @classmethod
    def get_builder_from_protocol(
        cls, pw_code, xs_code, structure, pw_protocol=None, overrides=None, xs_protocol=None, **kwargs
    ):
        """Return a builder prepopulated with inputs selected according to the chosen protocol.

        :param pw_code: the ``Code`` instance configured for the ``quantumespresso.pw``
            plugin.
        :param xs_code: the ``Code`` instance configured for the
            ``quantumespresso.xspectra`` plugin.
        :param structure: the ``StructureData`` instance to use.
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
            left=PwBaseWorkChain.get_protocol_inputs(protocol=pw_protocol,),
            right=cls.get_protocol_inputs(protocol=xs_protocol, overrides=overrides)
        )

        args = (pw_code, structure, pw_protocol)
        scf = PwBaseWorkChain.get_builder_from_protocol(*args, overrides=inputs.get('scf', None), **kwargs)

        scf.pop('clean_workdir', None)
        scf['pw'].pop('structure', None)
        # pylint: disable=no-member
        builder = cls.get_builder()
        builder.scf = scf
        builder.xs_prod.code = xs_code
        builder.xs_prod.parameters = orm.Dict(inputs.get('xs_prod', {}).get('parameters'))
        builder.xs_prod.metadata = inputs.get('xs_prod', {}).get('metadata')
        builder.xs_plot.code = xs_code
        builder.xs_plot.parameters = orm.Dict(inputs.get('xs_plot', {}).get('parameters'))
        builder.xs_plot.metadata = inputs.get('xs_prod', {}).get('metadata')
        builder.structure = structure
        builder.clean_workdir = orm.Bool(inputs['clean_workdir'])
        builder.collect_powder = orm.Bool(inputs['collect_powder'])
        builder.abs_atom_marker = orm.Str(inputs['abs_atom_marker'])
        # pylint: enable=no-member
        return builder

    def setup(self):
        """Initialize context variables that are used during the logical flow of the workchain."""

        self.ctx.serial_clean = 'serial_clean' in self.inputs and self.inputs.serial_clean.value
        self.ctx.all_lanczos_computed = False
        self.ctx.lanczos_to_restart = []
        self.ctx.finished_lanczos = []

    def should_repeat_xs_prod(self):
        """Return whether the Lanczos production step is finished or not."""

        return not self.ctx.all_lanczos_computed

    def run_scf(self):
        """Run an SCF calculation as a first step."""

        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, 'scf'))

        inputs.metadata.call_link_label = 'scf'
        inputs.pw.structure = self.inputs.structure
        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)

        future = self.submit(PwBaseWorkChain, **inputs)
        self.report(f'launching SCF PwBaseWorkChain<{future.pk}>')

        return ToContext(scf_workchain=future)

    def inspect_scf(self):
        """Verify that the PwBaseWorkChain finished successfully."""

        workchain = self.ctx.scf_workchain
        if not workchain.is_finished_ok:
            self.report(f'SCF PwBaseWorkChain failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF

        self.ctx.scf_parent_folder = workchain.outputs.remote_folder
        self.ctx.scf_kpoint_mesh = workchain.called[0].get_incoming().get_node_by_label('kpoints')

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
        code_label = 'upf2plotcore.sh@localhost'

        pw_inputs = self.exposed_inputs(PwBaseWorkChain, 'scf')
        absorbing_species = self.inputs.abs_atom_marker.value

        pseudo_dict = pw_inputs['pw']['pseudos']
        upf = pseudo_dict[absorbing_species]

        shell_inputs = {}

        # This part would require a new profile or attempt to run the shell script on a
        # remote machine in order to test, so we will leave it aside for now.
        # try:
        #     plotcore_code = orm.load_code(code_label)
        # except exceptions.NotExistent:
        #     self.report('No upf2plotcore.sh@localhost code found, creating one now')
        #     plotcore_code = orm.Code(
        #         label='upf2plotcore.sh',
        #         remote_computer_exec=('localhost', 'upf2plotcore.sh'),
        #         input_plugin_name='core.shell'
        #     ).store()

        # For now, the code node will be the one that we already have, but in the future we
        # will need to re-use the code in the launch_shell_job method in order to dynamically
        # produce the required code node in the same way.
        shell_inputs['code'] = orm.load_code(code_label)
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

        if 'core_wfc_data' in self.inputs:
            core_wfc_data = self.inputs.core_wfc_data
        else:
            core_wfc_data = self.ctx.upf2plotcore_node.get_outgoing().get_node_by_label('stdout')

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
        xspectra_calc_labels = []
        calc_number = 0
        if len(self.ctx.lanczos_to_restart) == 0: # No restarts have been ordered yet
            for vector in eps_vectors:
                xspectra_inputs = AttributeDict(self.exposed_inputs(XspectraCalculation, 'xs_prod'))
                xspectra_parameters = xspectra_inputs.parameters.get_dict()
                xspectra_parameters['INPUT_XSPECTRA']['xiabs'] = xiabs
                max_walltime_seconds = self.inputs.xs_prod.metadata.options.max_wallclock_seconds
                xspectra_parameters['INPUT_XSPECTRA']['time_limit'] = max_walltime_seconds * 0.9

                xspectra_inputs.parent_folder = self.ctx.scf_parent_folder
                xspectra_inputs.kpoints = self.ctx.scf_kpoint_mesh
                xspectra_inputs.core_wfc_data = core_wfc_data
                label = f'xas_prod_{calc_number}'
                xspectra_inputs.metadata.call_link_label = f'{label}_iter_1'
                xspectra_inputs.metadata.label = f'{label}_iter_1'
                xspectra_calc_labels.append(label)

                for index in [0, 1, 2]:
                    xspectra_parameters['INPUT_XSPECTRA'][f'xepsilon({index + 1})'] = vector[index]
                xspectra_inputs.parameters = orm.Dict(dict=xspectra_parameters)

                future_xspectra = self.submit(XspectraCalculation, **xspectra_inputs)
                self.report(
                    f'launching XspectraCalculation<{future_xspectra.pk}> for epsilon vector {vector}'
                    ' (Lanczos production) (iteration #1)'
                )
                xspectra_prod_calcs[label] = future_xspectra
                calc_number += 1

        else: # Some calculations need restarting, so we process these instead
            for calculation in self.ctx.lanczos_to_restart:
                xspectra_inputs = AttributeDict(self.exposed_inputs(XspectraCalculation, 'xs_prod'))
                xspectra_parameters = calculation.get_incoming().get_node_by_label('parameters').get_dict()
                xspectra_parameters['INPUT_XSPECTRA']['restart_mode'] = 'restart'
                vector = [
                    xspectra_parameters['INPUT_XSPECTRA']['xepsilon(1)'],
                    xspectra_parameters['INPUT_XSPECTRA']['xepsilon(2)'],
                    xspectra_parameters['INPUT_XSPECTRA']['xepsilon(3)']
                ]

                xspectra_inputs.parent_folder = calculation.get_outgoing().get_node_by_label('remote_folder')
                xspectra_inputs.kpoints = self.ctx.scf_kpoint_mesh
                xspectra_inputs.core_wfc_data = core_wfc_data
                parent_label_pieces = calculation.label.split('_')
                iteration = int(parent_label_pieces[-1])
                iteration += 1
                label = f'{parent_label_pieces[0]}_{parent_label_pieces[1]}_{parent_label_pieces[2]}'
                xspectra_inputs.metadata.call_link_label = f'{label}_iter_{iteration}'
                xspectra_inputs.metadata.label = f'{label}_iter_{iteration}'
                xspectra_calc_labels.append(label)

                xspectra_inputs.parameters = orm.Dict(dict=xspectra_parameters)

                future_xspectra = self.submit(XspectraCalculation, **xspectra_inputs)
                self.report(
                    f'launching XspectraCalculation<{future_xspectra.pk}> for epsilon vector {vector}'
                    ' (Lanczos production) (iteration #{iteration})'
                )
                xspectra_prod_calcs[f'{label}'] = future_xspectra

        self.ctx.xspectra_calc_labels = xspectra_calc_labels
        return ToContext(**xspectra_prod_calcs)

    def inspect_all_xspectra_prod(self):
        """Verify that the `XspectraCalculation` Lanczos production sub-processes finished successfully."""

        calculations = []
        labels = self.ctx.xspectra_calc_labels
        for label in labels:
            calculation = self.ctx[label]
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

        if 'core_wfc_data' in self.inputs:
            core_wfc_data = self.inputs.core_wfc_data
        else:
            core_wfc_data = self.ctx.upf2plotcore_node.get_outgoing().get_node_by_label('stdout')

        xspectra_plot_calcs = {}
        for calculation in finished_calculations:
            xspectra_inputs = AttributeDict(self.exposed_inputs(XspectraCalculation, 'xs_plot'))
            # The epsilon or k vectors are not needed in the case of a replot, however they
            # will be needed by the Parser at the end
            xspectra_parameters = xspectra_inputs.parameters.get_dict()

            parent_xspectra = calculation
            parent_output_dict = parent_xspectra.get_outgoing().get_node_by_label('output_parameters').get_dict()
            eps_vector = parent_output_dict['epsilon_vector']
            xspectra_parameters['INPUT_XSPECTRA']['xepsilon(1)'] = eps_vector[0]
            xspectra_parameters['INPUT_XSPECTRA']['xepsilon(2)'] = eps_vector[1]
            xspectra_parameters['INPUT_XSPECTRA']['xepsilon(3)'] = eps_vector[2]
            parent_label_pieces = calculation.label.split('_')
            label = f'{parent_label_pieces[0]}_{parent_label_pieces[1]}_{parent_label_pieces[2]}'
            xspectra_inputs.parent_folder = parent_xspectra.get_outgoing().get_node_by_label('remote_folder')
            xspectra_inputs.kpoints = self.ctx.scf_kpoint_mesh
            xspectra_inputs.core_wfc_data = core_wfc_data
            xspectra_inputs.metadata.label = label.replace('prod', 'plot')
            xspectra_inputs.metadata.call_link_label = label.replace('prod', 'plot')

            xspectra_inputs.parameters = orm.Dict(dict=xspectra_parameters)

            future_xspectra = self.submit(XspectraCalculation, **xspectra_inputs)
            self.report(f'launching XspectraCalculation<{future_xspectra.pk}> for epsilon vector {eps_vector} (Replot)')
            xspectra_plot_calcs[f'{label}_plot'] = future_xspectra

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

        This will collect the SCF and XSpectra output parameters,as well as the
        'Powder' spectrum (if requested)
        """

        xspectra_prod_calcs = self.ctx.finished_lanczos
        xspectra_plot_calcs = self.ctx.finished_replots

        should_collect_powder = False
        eps_powder_vectors = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
        for calc in xspectra_prod_calcs:
            output_parameters = calc.get_outgoing().get_node_by_label('output_parameters').get_dict()
            eps_vectors = output_parameters['epsilon_vector']
            if eps_vectors in eps_powder_vectors:
                should_collect_powder = True

        eps_basis_calcs = {}
        if self.inputs.collect_powder and should_collect_powder:
            for plot_calc in xspectra_plot_calcs:
                parent_folder = plot_calc.get_incoming().get_node_by_label('parent_folder')
                parent_calc = parent_folder.creator
                parent_output_params = parent_calc.get_outgoing().get_node_by_label('output_parameters').get_dict()
                parent_vector = parent_output_params['epsilon_vector']
                if parent_vector in eps_powder_vectors:
                    if parent_vector == [1., 0., 0.]:
                        eps_basis_calcs['eps_100'] = plot_calc.pk
                    if parent_vector == [0., 1., 0.]:
                        eps_basis_calcs['eps_010'] = plot_calc.pk
                    if parent_vector == [0., 0., 1.]:
                        eps_basis_calcs['eps_001'] = plot_calc.pk

            eps_basis_calcs['metadata'] = {'call_link_label': 'compile_powder'}
            powder_spectrum = get_powder_spectrum(**eps_basis_calcs)
            self.out('powder_spectrum', powder_spectrum)
        elif self.inputs.collect_powder and not should_collect_powder:
            # This should be upgraded to a more obvious warning, but I'm not sure yet how
            # to do this in practice.
            self.report(
                'Warning: A powder spectrum was requested, but none of the epsilon vectors '
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
            output_params = node.get_outgoing().get_node_by_label('output_parameters')
            xspectra_prod_params[key] = output_params
        self.out('output_parameters_xspectra', xspectra_prod_params)

        all_xspectra_plot_pks = {}
        for calc in xspectra_plot_calcs:
            all_xspectra_plot_pks[calc.label] = calc.pk

        all_xspectra_plot_pks['metadata'] = {'call_link_label': 'compile_all_spectra'}
        output_spectra = get_all_spectra(**all_xspectra_plot_pks)

        self.out('output_spectra', output_spectra)
        if self.inputs.clean_workdir.value is True:
            self.report('workchain succesfully completed, cleaning remote folders')
        elif self.inputs.clean_workdir.value is False:
            self.report('workchain succesfully completed, remote folders will be kept')

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
