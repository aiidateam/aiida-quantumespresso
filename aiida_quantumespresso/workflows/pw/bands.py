# -*- coding: utf-8 -*-
from aiida.common.extendeddicts import AttributeDict
from aiida.orm.calculation import JobCalculation
from aiida.orm.data.base import Str, Bool
from aiida.orm.data.parameter import ParameterData
from aiida.orm.data.structure import StructureData
from aiida.orm.data.array.bands import BandsData
from aiida.orm.group import Group
from aiida.orm.utils import WorkflowFactory
from aiida.work.workchain import WorkChain, ToContext, if_
from aiida.work.workfunctions import workfunction
from aiida_quantumespresso.utils.mapping import prepare_process_inputs


PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')
PwRelaxWorkChain = WorkflowFactory('quantumespresso.pw.relax')


class PwBandsWorkChain(WorkChain):
    """Workchain to compute a band structure for a given structure using Quantum ESPRESSO pw.x"""

    @classmethod
    def define(cls, spec):
        super(PwBandsWorkChain, cls).define(spec)
        spec.expose_inputs(PwRelaxWorkChain, namespace='relax', exclude=('structure', 'clean_workdir'))
        spec.expose_inputs(PwBaseWorkChain, namespace='scf', exclude=('structure', 'clean_workdir'))
        spec.expose_inputs(PwBaseWorkChain, namespace='bands', exclude=('structure', 'clean_workdir'))
        spec.input('structure', valid_type=StructureData)
        spec.input('clean_workdir', valid_type=Bool, default=Bool(False))
        spec.input('group', valid_type=Str, required=False)
        # These options are used to control the WorkChain behaviour. For instance,
        # to increase the number of empty bands to compute. They are checked in validate_inputs
        # (where you can check more in detail which ones are accepted).
        spec.input('workchain_options',  valid_type=ParameterData, required=False)
        spec.outline(
            cls.setup,
            cls.validate_inputs,
            if_(cls.should_do_relax)(
                cls.run_relax,
                cls.inspect_relax,
            ),
            cls.run_seekpath,
            cls.run_scf,
            cls.inspect_scf,
            cls.run_bands,
            cls.inspect_bands,
            cls.results,
        )
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_RELAX',
            message='the PwRelaxWorkChain sub process failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_SCF',
            message='the scf PwBasexWorkChain sub process failed')
        spec.exit_code(403, 'ERROR_SUB_PROCESS_FAILED_BANDS',
            message='the bands PwBasexWorkChain sub process failed')
        spec.exit_code(501, 'ERROR_INVALID_WORKCHAIN_OPTIONS_INPUTS',
            message='invalid inputs provided in workchain_options')
        spec.output('primitive_structure', valid_type=StructureData)
        spec.output('seekpath_parameters', valid_type=ParameterData)
        spec.output('scf_parameters', valid_type=ParameterData)
        spec.output('band_parameters', valid_type=ParameterData)
        spec.output('band_structure', valid_type=BandsData)

    def setup(self):
        """Define the current structure in the context to be the input structure."""
        self.ctx.current_structure = self.inputs.structure

    def validate_inputs(self):
        """Validate the additional inputs to this workchain. Make sure """
        try:
            wc_options = self.inputs['workchain_options'].get_dict()
        except KeyError:
            wc_options = {}

        # This value sets the number of bands for the band structure calculations
        # defined as a multiplicative factore wrt to the number of occupied states
        # e.g. 1.2 means 20% more or at least 4 more as in QE.
        # If not specified, set it to None. It will be then set to a proper
        # default inside run_bands
        self.ctx.num_bands_factor = wc_options.pop('num_bands_factor', None)

        if wc_options:
            #Checking that all options given in input are recognised
            self.report('Unknown variable passed inside workchain_options: {}'.format(
                              ', '.join(wc_options.keys())
                              ))
            return self.ERROR_INVALID_WORKCHAIN_OPTIONS_INPUTS

    def should_do_relax(self):
        """If the 'relax' input namespace was specified, we relax the input structure."""
        return 'relax' in self.inputs

    def run_relax(self):
        """Run the PwRelaxWorkChain to run a relax PwCalculation."""
        inputs = AttributeDict(self.exposed_inputs(PwRelaxWorkChain, namespace='relax'))
        inputs.structure = self.ctx.current_structure

        running = self.submit(PwRelaxWorkChain, **inputs)

        self.report('launching PwRelaxWorkChain<{}>'.format(running.pk))

        return ToContext(workchain_relax=running)

    def inspect_relax(self):
        """Verify that the PwRelaxWorkChain finished successfully."""
        workchain = self.ctx.workchain_relax

        if not workchain.is_finished_ok:
            self.report('PwRelaxWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_RELAX
        else:
            self.ctx.current_structure = workchain.out.output_structure

    def run_seekpath(self):
        """
        Run the relaxed structure through SeeKPath to get the new primitive structure, just in case
        the symmetry of the cell changed in the cell relaxation step
        """
        if 'kpoints_distance' in self.inputs.bands:
            seekpath_parameters = ParameterData(dict={
                'reference_distance': self.inputs.bands.kpoints_distance.value
            })
        else:
            seekpath_parameters = ParameterData(dict={})

        result = seekpath_structure_analysis(self.ctx.current_structure, seekpath_parameters)
        self.ctx.current_structure = result['primitive_structure']
        self.ctx.kpoints_path = result['explicit_kpoints']

        self.out('primitive_structure', result['primitive_structure'])
        self.out('seekpath_parameters', result['parameters'])

    def run_scf(self):
        """Run the PwBaseWorkChain in scf mode on the primitive cell of (optionally relaxed) input structure"""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='scf'))
        inputs.structure = self.ctx.current_structure
        inputs.parameters = inputs.parameters.get_dict()
        inputs.parameters.setdefault('CONTROL', {})
        inputs.parameters['CONTROL']['calculation'] = 'scf'
        
        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)

        running = self.submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}> in {} mode'.format(running.pk, 'scf'))

        return ToContext(workchain_scf=running)

    def inspect_scf(self):
        """Verify that the PwBaseWorkChain for the scf run finished successfully."""
        workchain = self.ctx.workchain_scf

        if not workchain.is_finished_ok:
            self.report('scf PwBaseWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SCF
        else:
            self.ctx.current_folder = workchain.out.remote_folder

    def run_bands(self):
        """Run the PwBaseWorkChain in bands mode along the path of high-symmetry determined by seekpath."""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain, namespace='bands'))
        inputs.parameters = inputs.parameters.get_dict()
        inputs.parameters.setdefault('CONTROL', {})
        inputs.parameters['CONTROL']['restart_mode'] = 'restart'
        inputs.parameters['CONTROL']['calculation'] = 'bands'
        # c-g seems better for bands as a default
        inputs.parameters['ELECTRONS']['diagonalization'] = 'cg'
        # We need full accuracy as we are doing bands
        inputs.parameters['ELECTRONS']['diago_full_acc'] = True

        ## START - Manage the number of additional bands ##
        # Get info from SCF on number of electrons and number of spin components
        scf_out_dict = self.ctx.workchain_scf.out.output_parameters.get_dict()
        num_elec = int(scf_out_dict['number_of_electrons'])
        num_spin = int(scf_out_dict['number_of_spin_components'])

        # This gives the same results also with noncollinear calcs
        # e.g. with 8 electrons, no spinors and a factor of 1.5 you compute 6 bands
        #      with spinors you compute 12 bands (still 50% more than the occupied states)
        # As in QE, we add at least 4 (8 with spinors) additional bands
        # If no num_bands_factor was specified (and the validation above set it to None),
        # use 1.2 as the default value (that is the default of QE as well: 20% more bands),
        # but do it also in the case of insulators.
        num_bands_factor = self.ctx.num_bands_factor
        if num_bands_factor is None:
            num_bands_factor = 1.2

        inputs.parameters['SYSTEM']['nbnd'] = max(
            int(0.5 * num_elec * num_spin * num_bands_factor),
            int(num_elec * num_spin * 0.5) + 4 * num_spin)
        ## END - Manage the number of additional bands ##

        if 'kpoints' not in self.inputs.bands:
            inputs.kpoints = self.ctx.kpoints_path

        inputs.structure = self.ctx.current_structure
        inputs.parent_folder = self.ctx.current_folder

        inputs = prepare_process_inputs(PwBaseWorkChain, inputs)
        running = self.submit(PwBaseWorkChain, **inputs)

        self.report('launching PwBaseWorkChain<{}> in {} mode'.format(running.pk, 'bands'))

        return ToContext(workchain_bands=running)

    def inspect_bands(self):
        """Verify that the PwBaseWorkChain for the bands run finished successfully."""
        workchain = self.ctx.workchain_bands

        if not workchain.is_finished_ok:
            self.report('bands PwBaseWorkChain failed with exit status {}'.format(workchain.exit_status))
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_BANDS

    def results(self):
        """Attach the desired output nodes directly as outputs of the workchain."""
        self.report('workchain succesfully completed')
        self.out('scf_parameters', self.ctx.workchain_scf.out.output_parameters)
        self.out('band_parameters', self.ctx.workchain_bands.out.output_parameters)
        self.out('band_structure', self.ctx.workchain_bands.out.output_band)

        if 'group' in self.inputs:
            output_band = self.ctx.workchain_bands.out.output_band
            group, _ = Group.get_or_create(name=self.inputs.group.value)
            group.add_nodes(output_band)
            self.report("storing the output_band<{}> in the group '{}'"
                .format(output_band.pk, self.inputs.group.value))

    def on_terminated(self):
        """
        If the clean_workdir input was set to True, recursively collect all called Calculations by
        ourselves and our called descendants, and clean the remote folder for the JobCalculation instances
        """
        super(PwBandsWorkChain, self).on_terminated()

        if self.inputs.clean_workdir.value is False:
            self.report('remote folders will not be cleaned')
            return

        cleaned_calcs = []

        for called_descendant in self.calc.called_descendants:
            if isinstance(called_descendant, JobCalculation):
                try:
                    called_descendant.out.remote_folder._clean()
                    cleaned_calcs.append(called_descendant.pk)
                except (IOError, OSError, KeyError):
                    pass

        if cleaned_calcs:
            self.report('cleaned remote folders of calculations: {}'.format(' '.join(map(str, cleaned_calcs))))


@workfunction
def seekpath_structure_analysis(structure, parameters):
    """
    This workfunction will take a structure and pass it through SeeKpath to get the
    primitive cell and the path of high symmetry k-points through its Brillouin zone.
    Note that the returned primitive cell may differ from the original structure in
    which case the k-points are only congruent with the primitive cell.
    """
    from aiida.tools import get_explicit_kpoints_path
    return get_explicit_kpoints_path(structure, **parameters.get_dict())
