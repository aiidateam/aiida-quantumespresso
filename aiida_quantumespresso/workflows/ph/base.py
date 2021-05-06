# -*- coding: utf-8 -*-
"""Workchain to run a Quantum ESPRESSO ph.x calculation with automated error handling and restarts."""
from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import while_, BaseRestartWorkChain, process_handler, ProcessHandlerReport, calcfunction
from aiida.plugins import CalculationFactory

PhCalculation = CalculationFactory('quantumespresso.ph')
PwCalculation = CalculationFactory('quantumespresso.pw')


@calcfunction
def merge_outputs(**kwargs):
    #take output from last run
    #from other runs merge the parsed q-points &  found irreps,
    #                add times num_q_found
    
    final_index = max(kwargs.keys())
    merged = kwargs[final_index].get_dict()
    merged['number_of_irr_representations_for_each_q'] = [None for i in range(merged['number_of_qpoints'])]
    for index, child_dict in kwargs.items():
        which_q=[]
        if index == final_index:
            for k in child_dict.keys():
                if 'q_point_' in k and 'mode_symmetry' in child_dict[k].keys():
                    which_q.append(int(k.split('_')[-1]))
        else:
            for k in child_dict.keys():
                if k == 'num_q_found' or k == 'wall_time_seconds':
                    merged[k] += child_dict[k]
                elif 'q_point_' in k and 'mode_symmetry' in child_dict[k].keys():
                    #this qpoint has been parsed - keep
                    merged[k] = child_dict[k]
                    which_q.append(int(k.split('_')[-1]))

        which_q.sort()
        for i,qind in enumerate(which_q):
            merged['number_of_irr_representations_for_each_q'][qind-1] = child_dict['number_of_irr_representations_for_each_q'][i]
        
    merged['wall_time'] = str(merged['wall_time_seconds']//60)+'m'+str(merged['wall_time_seconds']-merged['wall_time_seconds']//60*60)+'s'

    merged = orm.Dict(dict=merged)
    return merged
        


class PhBaseWorkChain(BaseRestartWorkChain):
    """Workchain to run a Quantum ESPRESSO ph.x calculation with automated error handling and restarts."""

    _process_class = PhCalculation

    defaults = AttributeDict({
        'delta_factor_max_seconds': 0.95,
        'delta_factor_alpha_mix': 0.90,
        'alpha_mix': 0.70,
    })

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.expose_inputs(PhCalculation, namespace='ph')
        spec.input('only_initialization', valid_type=orm.Bool, default=lambda: orm.Bool(False))
        spec.outline(
            cls.setup,
            cls.validate_parameters,
            cls.validate_resources,
            while_(cls.should_run_process)(
                cls.prepare_process,
                cls.run_process,
                cls.inspect_process,
            ),
            cls.results,
            cls.create_merged_output,
        )
        #spec.expose_outputs(PhCalculation, exclude=('output_parameters',))


        spec.expose_outputs(PwCalculation, exclude=('retrieved_folder',))
        spec.output('merged_output_parameters', valid_type=orm.Dict, required=False)
        spec.exit_code(204, 'ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED',
            message='The `metadata.options` did not specify both `resources.num_machines` and `max_wallclock_seconds`.')
        spec.exit_code(300, 'ERROR_UNRECOVERABLE_FAILURE',
            message='The calculation failed with an unrecoverable error.')
        # yapf: enable

    def setup(self):
        """Call the `setup` of the `BaseRestartWorkChain` and then create the inputs dictionary in `self.ctx.inputs`.

        This `self.ctx.inputs` dictionary will be used by the `BaseRestartWorkChain` to submit the calculations in the
        internal loop.
        """
        super().setup()
        self.ctx.restart_calc = None
        self.ctx.inputs = AttributeDict(self.exposed_inputs(PhCalculation, 'ph'))
        
    def validate_parameters(self):
        """Validate inputs that might depend on each other and cannot be validated by the spec."""
        self.ctx.inputs.parameters = self.ctx.inputs.parameters.get_dict()
        self.ctx.inputs.settings = self.ctx.inputs.settings.get_dict() if 'settings' in self.ctx.inputs else {}

        self.ctx.inputs.parameters.setdefault('INPUTPH', {})
        self.ctx.inputs.parameters['INPUTPH']['recover'] = 'parent_folder' in self.ctx.inputs

        if self.inputs.only_initialization.value:
            self.ctx.inputs.settings['ONLY_INITIALIZATION'] = True

    def validate_resources(self):
        """Validate the inputs related to the resources.

        The `metadata.options` should at least contain the options `resources` and `max_wallclock_seconds`, where
        `resources` should define the `num_machines`.
        """
        num_machines = self.ctx.inputs.metadata.options.get('resources', {}).get('num_machines', None)
        max_wallclock_seconds = self.ctx.inputs.metadata.options.get('max_wallclock_seconds', None)

        if num_machines is None or max_wallclock_seconds is None:
            return self.exit_codes.ERROR_INVALID_INPUT_RESOURCES_UNDERSPECIFIED

        self.set_max_seconds(max_wallclock_seconds)

    def set_max_seconds(self, max_wallclock_seconds):
        """Set the `max_seconds` to a fraction of `max_wallclock_seconds` option to prevent out-of-walltime problems.

        :param max_wallclock_seconds: the maximum wallclock time that will be set in the scheduler settings.
        """
        max_seconds_factor = self.defaults.delta_factor_max_seconds
        max_seconds = min(max_wallclock_seconds - 60, max_wallclock_seconds * max_seconds_factor)
        self.ctx.inputs.parameters['INPUTPH']['max_seconds'] = max_seconds

    def prepare_process(self):
        """Prepare the inputs for the next calculation.

        If a `restart_calc` has been set in the context, its `remote_folder` will be used as the `parent_folder` input
        for the next calculation and the `restart_mode` is set to `restart`.
        """
        if self.ctx.restart_calc:
            self.ctx.inputs.parameters['INPUTPH']['recover'] = True
            self.ctx.inputs.parent_folder = self.ctx.restart_calc.outputs.remote_folder


    def create_merged_output(self):
        outputs = {}
        for index, child in enumerate(self.ctx.children):
            outputs['child_'+str(index)] = child.outputs.output_parameters

        #check that runs have the same number of qpoints & all were parsed
        num_q = None
        q_found = 0
        for index, child_dict in outputs.items():
            q_found += child_dict['num_q_found']
            if num_q:
                if num_q == child_dict['number_of_qpoints']:
                    pass
                else:
                    #need error
                    self.report('All PhCalculations do not have the same number of q-points')
            else:
                num_q = child_dict['number_of_qpoints']
                    
        if q_found == num_q:
            self.report('Merging {} q-points'.format(q_found))
            self.out('merged_output_parameters', merge_outputs(**outputs))
        else: #need error
            self.report('Only {} of {} q-points were parsed'.format(q_found,num_q))
            
        
 
                                  
    def report_error_handled(self, calculation, action):
        """Report an action taken for a calculation that has failed.

        This should be called in a registered error handler if its condition is met and an action was taken.

        :param calculation: the failed calculation node
        :param action: a string message with the action taken
        """
        arguments = [calculation.process_label, calculation.pk, calculation.exit_status, calculation.exit_message]
        self.report('{}<{}> failed with exit status {}: {}'.format(*arguments))
        self.report(f'Action taken: {action}')

    @process_handler(priority=600)
    def handle_unrecoverable_failure(self, node):
        """Handle calculations with an exit status below 400 which are unrecoverable, so abort the work chain."""
        if node.is_failed and node.exit_status < 400:
            self.report_error_handled(node, 'unrecoverable error, aborting...')
            return ProcessHandlerReport(True, self.exit_codes.ERROR_UNRECOVERABLE_FAILURE)

    @process_handler(priority=580, exit_codes=PhCalculation.exit_codes.ERROR_OUT_OF_WALLTIME)
    def handle_out_of_walltime(self, node):
        """Handle `ERROR_OUT_OF_WALLTIME` exit code: calculation shut down neatly and we can simply restart."""
        self.ctx.restart_calc = node
        self.report_error_handled(node, 'simply restart from the last calculation')
        return ProcessHandlerReport(True)

    @process_handler(priority=410, exit_codes=PhCalculation.exit_codes.ERROR_CONVERGENCE_NOT_REACHED)
    def handle_convergence_not_achieved(self, node):
        """Handle `ERROR_CONVERGENCE_NOT_REACHED` exit code: decrease the mixing beta and restart from scratch."""
        factor = self.defaults.delta_factor_alpha_mix
        alpha_mix = self.ctx.inputs.parameters.get('INPUTPH', {}).get('alpha_mix(1)', self.defaults.alpha_mix)
        alpha_mix_new = alpha_mix * factor

        self.ctx.restart_calc = node
        self.ctx.inputs.parameters.setdefault('INPUTPH', {})['alpha_mix(1)'] = alpha_mix_new

        action = f'reduced alpha_mix from {alpha_mix} to {alpha_mix_new} and restarting'
        self.report_error_handled(node, action)
        return ProcessHandlerReport(True)
