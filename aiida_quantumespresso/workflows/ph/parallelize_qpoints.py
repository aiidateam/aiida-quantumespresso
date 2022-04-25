# -*- coding: utf-8 -*-
"""Workchain to perform a ph.x calculation with parallelization over q-points."""
from aiida import orm
from aiida.common.extendeddicts import AttributeDict
from aiida.engine import ToContext, WorkChain, append_, calcfunction
from aiida.plugins import CalculationFactory, WorkflowFactory

from aiida_quantumespresso.utils.mapping import prepare_process_inputs
from aiida_quantumespresso.utils.resources import get_default_options

PhCalculation = CalculationFactory('quantumespresso.ph')
PwCalculation = CalculationFactory('quantumespresso.pw')
PhBaseWorkChain = WorkflowFactory('quantumespresso.ph.base')


class PhParallelizeQpointsWorkChain(WorkChain):
    """Workchain to launch a Quantum Espresso phonon ph.x calculation.

    This workchain differs from the `PhBaseWorkChain` in that the computation is parallelized over the q-points. For
    each individual q-point a separate `PhBaseWorkChain` is run, distributed over the provided resources.
    """

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        super().define(spec)
        spec.input('code', valid_type=orm.Code)
        spec.input('parent_folder', valid_type=orm.RemoteData)
        spec.input('qpoints', valid_type=orm.KpointsData)
        spec.input('parameters', valid_type=orm.Dict, required=False)
        spec.input('settings', valid_type=orm.Dict, required=False)
        spec.input('options', valid_type=orm.Dict, required=False)
        spec.input('clean_workdir', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.input('compute_epsil', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.input('alpha_mix', valid_type=orm.Float, default=lambda: orm.Float(0.7))
        spec.input('max_iterations', valid_type=orm.Int, default=lambda: orm.Int(10))
        spec.outline(
            cls.validate_inputs,
            cls.run_ph_init,
            cls.run_ph_qgrid,
            cls.results,
            cls.run_clean,
        )
        spec.output('retrieved', valid_type=orm.FolderData)

    def validate_inputs(self):
        """Validate inputs that depend might depend on each other and cannot be validated by the spec.

        Also define dictionary `inputs` in the context, that will contain the inputs for the sub workchain that will be
        launched.
        """
        self.ctx.inputs = AttributeDict({
            'code': self.inputs.code,
            'qpoints': self.inputs.qpoints,
            'parent_folder': self.inputs.parent_folder,
            'clean_workdir': self.inputs.clean_workdir,
            'compute_epsil': self.inputs.compute_epsil,
            'alpha_mix': self.inputs.alpha_mix,
            'max_iterations': self.inputs.max_iterations,
        })

        if 'parameters' in self.inputs:
            self.ctx.inputs.parameters = self.inputs.parameters.get_dict()
        else:
            self.ctx.inputs.parameters = {}

        if 'INPUTPH' not in self.ctx.inputs.parameters:
            self.ctx.inputs.parameters['INPUTPH'] = {}

        if 'settings' in self.inputs:
            self.ctx.inputs.settings = self.inputs.settings

        if 'options' in self.inputs:
            self.ctx.inputs.options = self.inputs.options

    def run_ph_init(self):
        """Run a first dummy ph calculation that will exit straight after initialization.

        At that point it will have generated the q-point list, which we use to determine how to distribute these over
        the available computational resources.
        """
        # Toggle the only initialization flag and define default options for minimal resources
        self.ctx.inputs.only_initialization = orm.Bool(True)
        self.ctx.inputs.options = get_default_options()

        inputs = prepare_process_inputs(PhBaseWorkChain, self.ctx.inputs)
        running = self.submit(PhBaseWorkChain, **inputs)

        self.report(f'launching initialization PhBaseWorkChain<{running.pk}>')

        return ToContext(ph_init=running)

    def run_ph_qgrid(self):
        """Distribute the q-points and launch individual PhBaseWorkChains for each q-point."""
        try:
            parent_calc = self.inputs.parent_folder.get_inputs(node_type=PwCalculation)[0]
        except IndexError:
            self.abort_nowait(
                f'could not retrieve the parent PwCalculation from parent_folder<{self.inputs.parent_folder.pk}>'
            )

        try:
            # If the parent has an output structure, i.e. vc-relax, we need that one
            structure = parent_calc.out.output_structure
        except AttributeError:
            # Otherwise, we take the input structure of the parent calculation
            structure = parent_calc.inp.structure

        # Remove the initialization only from the inputs
        self.ctx.inputs.pop('only_initialization')

        self.report('Distributing q-points for parallelization')
        retrieved = self.ctx.ph_init.out.retrieved
        qpoints_dict = distribute_qpoints(retrieved=retrieved, structure=structure)
        self.report('Distributing q-points completed')

        # Store the original compute_epsil value, either set by user in parameters or default to input value
        compute_epsil = self.ctx.inputs.parameters['INPUTPH'].get('epsil', self.inputs.compute_epsil.value)

        for k in sorted(qpoints_dict.keys()):

            qpoint = qpoints_dict[k]

            if not all(_ == 0. for _ in qpoint.get_kpoints()[0]):
                # Can only compute the dielectric constant at Gamma
                self.ctx.inputs.parameters['INPUTPH']['epsil'] = False
            else:
                # Restore the parameter set by the user/default
                self.ctx.inputs.parameters['INPUTPH']['epsil'] = compute_epsil

            self.ctx.inputs.qpoints = qpoint

            inputs = prepare_process_inputs(PhBaseWorkChain, self.ctx.inputs)
            running = self.submit(PhBaseWorkChain, **inputs)

            self.report(f'launching PhBaseWorkChain<{running.pk}> for q-point<{qpoint.pk}>')
            self.to_context(workchains=append_(running))

    def results(self):
        """Collect all retrieved folders of the launched PhBaseWorkChains and merge them into a single FolderData."""
        retrieved_folders = {'0': self.ctx.ph_init.out.retrieved}

        for ph_workchain in self.ctx.workchains:
            qpoint_link = ph_workchain.inp.qpoints.get_inputs_dict().keys()[0]
            qpoint_index = str(int(qpoint_link.split('_')[1]) + 1)
            retrieved_folders[qpoint_index] = ph_workchain.out.retrieved

        merged_retrieved = recollect_qpoints(**retrieved_folders)

        self.out('retrieved', merged_retrieved)
        self.report('workchain completed successfully')

    def run_clean(self):
        """Clean remote folders of the initialization PhCalculation.

        The PhCalculations launched by the subworkchains will be handled by them.
        """
        if not self.inputs.clean_workdir.value:
            self.report('remote folders will not be cleaned')
            return

        try:
            calc = self.ctx.ph_init
            calc.remote_folder._clean()  # pylint: disable=protected-access
            self.report(f'cleaned remote folder of {calc.__class__.__name__}<{calc.pk}>')
        except Exception:  # pylint: disable=broad-except
            pass


@calcfunction
def distribute_qpoints(retrieved, structure):
    """Distribute qpoints mesh over individual k-points.

    :param retrieved: a FolderData object with the retrieved node of a PhCalculation
    :param structure: a StructureData object from the parent PwCalculation
    :return: KpointsData objects with link labels of form 'qpoint_N' where N is the q-point index
    """
    import numpy

    if not isinstance(structure, orm.StructureData):
        raise ValueError('the structure argument should be a StructureData object')

    if not isinstance(retrieved, orm.FolderData):
        raise ValueError('the retrieved argument should be a FolderData object')

    dynmat_prefix = PhCalculation()._OUTPUT_DYNAMICAL_MATRIX_PREFIX  # pylint: disable=protected-access
    dynmat_file = retrieved.get_abs_path(f'{dynmat_prefix}{0}')

    with open(dynmat_file, 'r', encoding='utf-8') as handle:
        lines = handle.readlines()

    try:
        _ = [float(i) for i in lines[0].split()]
    except (IndexError, ValueError) as exception:
        raise ValueError(f"File '{dynmat_file}' does not contain the list of q-points") from exception

    cell = structure.cell
    alat = numpy.linalg.norm(cell[0])
    fact = 2. * numpy.pi / alat

    # Read q-points, converting them from 2pi/a coordinates to inverse angstrom
    qpoint_coordinates = [[float(i) * fact for i in j.split()] for j in lines[2:]]

    qpoints = {}
    for index, qpoint_coordinate in enumerate(qpoint_coordinates):
        qpoint = orm.KpointsData()
        qpoint.set_cell(cell)
        qpoint.set_kpoints([qpoint_coordinate], cartesian=True)
        qpoints[f'qpoint_{index}'] = qpoint

    return qpoints


@calcfunction
def recollect_qpoints(**kwargs):
    """Collect dynamical matrix files into a single folder.

    For each dynamical matrix, a different number is put at the end of the file, obtained from the input link, which
    corresponds to its place in the list of q-points originally generated by distribute_qpoints.

    :param kwargs: keys are the string representation of the q-point index and the value is the corresponding retrieved
        folder object. A special case is the folder at key '0' which is the folder of the initialization calculation.
    :return: FolderData object containing the dynamic matrix files of the computed PhBaseWorkChains
    """
    dynmat_prefix = PhCalculation()._OUTPUT_DYNAMICAL_MATRIX_PREFIX  # pylint: disable=protected-access

    # Initialize the merged folder, by creating the subdirectory for the dynamical matrix files
    merged_folder = orm.FolderData()

    for key, retrieved_folder in kwargs.items():
        filepath_src = f'{dynmat_prefix}{key}'
        filepath_dst = f'{dynmat_prefix}{key}'

        # The dynamic matrix source file only has an index if it is from the initialization
        # calculation. Otherwise, the filename is just the prefix
        if key != '0':
            filepath_src = f'{dynmat_prefix}'

        with retrieved_folder.open(filepath_src, 'rb') as handle:
            merged_folder.put_object_from_filelike(handle, filepath_dst, mode='wb')

    return merged_folder
