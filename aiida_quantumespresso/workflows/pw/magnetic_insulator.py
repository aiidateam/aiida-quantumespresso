# -*- coding: utf-8 -*-
"""
Workchain to perform the double scf run for magnetic insulators
using Quantum ESPRESSO pw.x.
"""
from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import WorkChain, ToContext
from aiida.plugins import WorkflowFactory
from aiida.orm.nodes.data.array.bands import find_bandgap

from aiida_quantumespresso.utils.defaults.calculation import pw as qe_defaults
from aiida_quantumespresso.utils.validation.magnetization import set_tot_magnetization


PwBaseWorkChain = WorkflowFactory('quantumespresso.pw.base')


class MagneticInsulatoWorkChain(WorkChain):
    """
    Workchain to perform the double scf run for magnetic insulators
    using Quantum ESPRESSO pw.x.
    """
    
    defaults = AttributeDict({
        'qe': qe_defaults,
        'smearing_method': 'marzari-vanderbilt',
        'smearing_degauss': 0.01,
    })

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        
        spec.expose_inputs(PwBaseWorkChain)
        
        spec.outline(
            cls.run_scf_smearing,
            cls.inspect_scf_smearing,
            cls.run_scf_fixed,
            cls.inspect_scf_fixed,
            cls.results,
        )
        spec.exit_code(401, 'ERROR_SUB_PROCESS_FAILED_SMEARING',
            message='the scf PwBaseWorkChain with smearing occupations failed')
        spec.exit_code(402, 'ERROR_SUB_PROCESS_FAILED_FIXED',
            message='the final scf PwBaseWorkChain with fixed occupations failed')
        spec.expose_outputs(PwBaseWorkChain)
        # yapf: enable

    def run_scf_smearing(self):
        """Run a `PwBaseWorkChain` with smearing occupations."""
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain))
        parameters = inputs.pw.parameters.get_dict()
        
        parameters['CONTROL']['calculation'] = 'scf'
        parameters['SYSTEM']['occupations'] = 'smearing'
        parameters['SYSTEM']['smearing'] = parameters['SYSTEM'].get(
            'smearing', self.defaults.smearing_method
        )
        parameters['SYSTEM']['degauss'] = parameters['SYSTEM'].get(
            'degauss', self.defaults.smearing_degauss
        )
        
        inputs.pw.parameters = orm.Dict(dict=parameters)
        
        inputs.metadata.call_link_label = 'scf_smearing'
        running = self.submit(PwBaseWorkChain, **inputs)
        self.report(f'launching PwBaseWorkChain<{running.pk}> with smearing occupations')

        return ToContext(workchain_scf_smearing=running)

    def inspect_scf_smearing(self):
        """Inspect the result of the smearing occupations scf `PwBaseWorkChain`."""
        workchain = self.ctx.workchain_scf_smearing

        if not workchain.is_finished_ok:
            self.report(f'scf PwBaseWorkChain with smearing occupations failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_SMEARING
        
        bands = workchain.outputs.output_band
        parameters = workchain.outputs.output_parameters.get_dict()
        number_electrons = parameters['number_of_electrons']

        is_insulator, _ = find_bandgap(bands, number_electrons=number_electrons)

        if is_insulator:
            self.report('the ground state of the scf with smearing occupations '
            'has been recognized as insulating')
        else:
            self.report('WARNING - the ground state with smearing occupations '
            'has been recognized as metallic')
            
    def run_scf_fixed(self):
        """Run a `PwBaseWorkChain` with fixed occupations."""
        previous_workchain = self.ctx.workchain_scf_smearing
        previous_parameters = previous_workchain.outputs.output_parameters
        
        inputs = AttributeDict(self.exposed_inputs(PwBaseWorkChain))
        parameters = inputs.pw.parameters.get_dict()
        
        parameters['CONTROL']['calculation'] = 'scf'
        parameters['CONTROL']['restart_mode'] = 'from_scratch'
        parameters['SYSTEM']['occupations'] = 'fixed'
        
        parameters['SYSTEM'].pop('degauss', None)
        parameters['SYSTEM'].pop('smearing', None)
        parameters['SYSTEM'].pop('starting_magnetization', None)
        
        parameters['SYSTEM']['nbnd'] = previous_parameters.get_dict()['number_of_bands']
        
        if set_tot_magnetization( parameters, previous_parameters.get_dict()['total_magnetization'] ):
                return self.exit_codes.ERROR_NON_INTEGER_TOT_MAGNETIZATION.format(iteration=self.ctx.iteration)
        
        inputs.pw.parameters = orm.Dict(dict=parameters)
        inputs.pw.parent_folder = previous_workchain.outputs.remote_folder
        
        inputs.metadata.call_link_label = 'scf_fixed'
        running = self.submit(PwBaseWorkChain, **inputs)
        self.report(f'launching PwBaseWorkChain<{running.pk}> with fixed occupations')

        return ToContext(workchain_scf_fixed=running)

    def inspect_scf_fixed(self):
        """Inspect the result of the fixed occupations scf `PwBaseWorkChain`."""
        workchain = self.ctx.workchain_scf_fixed

        if not workchain.is_finished_ok:
            self.report(f'scf PwBaseWorkChain with fixed occupations failed with exit status {workchain.exit_status}')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED_FIXED

    def results(self):
        """Attach the output of the workchain with the fixed occupations to the outputs."""
        self.out_many(self.exposed_outputs(self.ctx.workchain_scf_fixed, PwBaseWorkChain))
