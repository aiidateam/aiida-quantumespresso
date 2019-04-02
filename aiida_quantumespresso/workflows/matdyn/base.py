# -*- coding: utf-8 -*-
from __future__ import absolute_import

from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import while_
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.common.workchain.base.restart import BaseRestartWorkChain
from aiida_quantumespresso.data.forceconstants import ForceconstantsData
from aiida_quantumespresso.utils.resources import get_default_options


MatdynCalculation = CalculationFactory('quantumespresso.matdyn')


class MatdynBaseWorkChain(BaseRestartWorkChain):
    """
    Base Workchain to launch a Quantum Espresso matdyn.x calculation and restart it until
    successfully finished or until the maximum number of restarts is exceeded
    """
    _calculation_class = MatdynCalculation

    @classmethod
    def define(cls, spec):
        super(MatdynBaseWorkChain, cls).define(spec)
        spec.input('code', valid_type=orm.Code)
        spec.input('kpoints', valid_type=orm.KpointsData)
        spec.input('parent_folder', valid_type=ForceconstantsData)
        spec.input('parameters', valid_type=orm.Dict, required=False)
        spec.input('settings', valid_type=orm.Dict, required=False)
        spec.input('options', valid_type=orm.Dict, required=False)
        spec.outline(
            cls.setup,
            cls.validate_inputs,
            while_(cls.should_run_calculation)(
                cls.run_calculation,
                cls.inspect_calculation,
            ),
            cls.results,
        )
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('output_phonon_bands', valid_type=orm.BandsData)

    def validate_inputs(self):
        """
        Validate inputs that depend might depend on each other and cannot be validated by the spec. Also define
        dictionary `inputs` in the context, that will contain the inputs for the calculation that will be launched
        in the `run_calculation` step.
        """
        self.ctx.inputs = AttributeDict({
            'code': self.inputs.code,
            'kpoints': self.inputs.kpoints,
            'parent_folder': self.inputs.parent_folder,
        })

        if 'parameters' in self.inputs:
            self.ctx.inputs.parameters = self.inputs.parameters.get_dict()
        else:
            self.ctx.inputs.parameters = {'INPUT': {}}

        if 'settings' in self.inputs:
            self.ctx.inputs.settings = self.inputs.settings.get_dict()
        else:
            self.ctx.inputs.settings = {}

        if 'options' in self.inputs:
            self.ctx.inputs.options = self.inputs.options.get_dict()
        else:
            self.ctx.inputs.options = get_default_options()
