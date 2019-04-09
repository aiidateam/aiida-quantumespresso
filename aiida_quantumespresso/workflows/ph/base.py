# -*- coding: utf-8 -*-
from __future__ import absolute_import

from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import while_
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.common.workchain.utils import ErrorHandlerReport, register_error_handler
from aiida_quantumespresso.common.workchain.base.restart import BaseRestartWorkChain
from aiida_quantumespresso.utils.resources import get_default_options


PhCalculation = CalculationFactory('quantumespresso.ph')
PwCalculation = CalculationFactory('quantumespresso.pw')


class PhBaseWorkChain(BaseRestartWorkChain):
    """
    Base Workchain to launch a Quantum Espresso phonon ph.x calculation and restart it until
    successfully converged or until the maximum number of restarts is exceeded
    """
    _calculation_class = PhCalculation

    defaults = AttributeDict({
        'delta_factor_max_seconds': 0.95,
        'alpha_mix': 0.70,
    })

    @classmethod
    def define(cls, spec):
        super(PhBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=orm.Code)
        spec.input('qpoints', valid_type=orm.KpointsData)
        spec.input('parent_folder', valid_type=orm.RemoteData)
        spec.input('parameters', valid_type=orm.Dict, required=False)
        spec.input('settings', valid_type=orm.Dict, required=False)
        spec.input('options', valid_type=orm.Dict, required=False)
        spec.input('only_initialization', valid_type=orm.Bool, default=orm.Bool(False))
        spec.outline(
            cls.setup,
            cls.validate_inputs,
            while_(cls.should_run_calculation)(
                cls.prepare_calculation,
                cls.run_calculation,
                cls.inspect_calculation,
            ),
            cls.results,
        )
        spec.exit_code(402, 'ERROR_CALCULATION_INVALID_INPUT_FILE',
            message='the calculation failed because it had an invalid input file')
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('remote_folder', valid_type=orm.RemoteData)
        spec.output('retrieved', valid_type=orm.FolderData)

    def validate_inputs(self):
        """
        Validate inputs that depend might depend on each other and cannot be validated by the spec. Also define
        dictionary `inputs` in the context, that will contain the inputs for the calculation that will be launched
        in the `run_calculation` step.
        """
        self.ctx.inputs = AttributeDict({
            'code': self.inputs.code,
            'qpoints': self.inputs.qpoints,
            'parent_folder': self.inputs.parent_folder,
        })

        if 'parameters' in self.inputs:
            self.ctx.inputs.parameters = self.inputs.parameters.get_dict()
        else:
            self.ctx.inputs.parameters = {}

        if 'INPUTPH' not in self.ctx.inputs.parameters:
            self.ctx.inputs.parameters['INPUTPH'] = {}

        if 'settings' in self.inputs:
            self.ctx.inputs.settings = self.inputs.settings.get_dict()
        else:
            self.ctx.inputs.settings = {}

        if 'options' in self.inputs:
            self.ctx.inputs.options = self.inputs.options.get_dict()
        else:
            self.ctx.inputs.options = get_default_options()

        if self.inputs.only_initialization.value:
            self.ctx.inputs.settings['ONLY_INITIALIZATION'] = True

    def prepare_calculation(self):
        """
        Prepare the inputs for the next calculation
        """
        if isinstance(self.ctx.restart_calc, PhCalculation):
            self.ctx.inputs.parameters['INPUTPH']['recover'] = True
            self.ctx.inputs.parent_folder = self.ctx.restart_calc.outputs.remote_folder

    def _prepare_process_inputs(self, process, inputs):
        """
        The 'max_seconds' setting in the 'INPUTPH' card of the parameters will be set to a fraction of the
        'max_wallclock_seconds' that will be given to the job via the 'options' dictionary. This will prevent the job
        from being prematurely terminated by the scheduler without getting the chance to exit cleanly.
        """
        max_wallclock_seconds = inputs.options['max_wallclock_seconds']
        max_seconds_factor = self.defaults.delta_factor_max_seconds
        max_seconds = max_wallclock_seconds * max_seconds_factor
        inputs.parameters['INPUTPH']['max_seconds'] = max_seconds

        return super(PhBaseWorkChain, self)._prepare_process_inputs(process, inputs)


@register_error_handler(PhBaseWorkChain, 400)
def _handle_fatal_error_read_namelists(self, calculation):
    """
    The calculation failed because it could not read the generated input file
    """
    if any(['reading inputph namelist' in w for w in calculation.res.warnings]):
        self.report('PhCalculation<{}> failed because of an invalid input file'.format(calculation.pk))
        return ErrorHandlerReport(True, True, self.exit_codes.ERROR_CALCULATION_INVALID_INPUT_FILE)


@register_error_handler(PhBaseWorkChain, 300)
def _handle_error_exceeded_maximum_walltime(self, calculation):
    """
    Calculation ended nominally but ran out of allotted wall time
    """
    if 'Maximum CPU time exceeded' in calculation.res.warnings:
        self.ctx.restart_calc = calculation
        self.report('PhCalculation<{}> terminated because maximum wall time was exceeded, restarting'
            .format(calculation.pk))
        return ErrorHandlerReport(True, True)


@register_error_handler(PhBaseWorkChain, 200)
def _handle_fatal_error_not_converged(self, calculation):
    """
    The calculation failed because it could not read the generated input file
    """
    if ('Phonon did not reach end of self consistency' in calculation.res.warnings):
        alpha_mix_old = calculation.inputs.parameters.get_dict()['INPUTPH'].get('alpha_mix(1)', self.defaults.alpha_mix)
        alpha_mix_new = 0.9 * alpha_mix_old
        self.ctx.inputs.parameters['INPUTPH']['alpha_mix(1)'] = alpha_mix_new
        self.ctx.restart_calc = calculation
        self.report('PhCalculation<{}> terminated without reaching convergence, '
            'setting alpha_mix to {} and restarting'.format(calculation.pk, alpha_mix_new))
        return ErrorHandlerReport(True, True)


@register_error_handler(PhBaseWorkChain, 100)
def _handle_error_premature_termination(self, calculation):
    """
    Calculation did not reach the end of execution, probably because it was killed by the scheduler
    for running out of allotted walltime
    """
    if 'QE ph run did not reach the end of the execution.' in calculation.res.parser_warnings:
        inputs = calculation.inputs.parameters.get_dict()
        settings = self.ctx.inputs.settings

        factor = self.defaults.delta_factor_max_seconds
        max_seconds = settings.get('max_seconds', inputs['INPUTPH']['max_seconds'])
        max_seconds_reduced = int(max_seconds * factor)
        self.ctx.inputs.parameters['INPUTPH']['max_seconds'] = max_seconds_reduced

        self.ctx.restart_calc = calculation
        self.report('PwCalculation<{}> was terminated prematurely, reducing "max_seconds" from {} to {}'
            .format(calculation.pk, max_seconds, max_seconds_reduced))
        return ErrorHandlerReport(True, True)
