# -*- coding: utf-8 -*-
import click
from aiida_quantumespresso.utils.click import command
from aiida_quantumespresso.utils.click import options

@command()
@options.code()
@options.parent_calc(callback_kwargs={'entry_point': 'quantumespresso.ph'})
@options.max_num_machines()
@options.max_wallclock_seconds()
def launch(
    code, parent_calc, max_num_machines, max_wallclock_seconds):
    """
    Run the Q2rBaseWorkChain for a previously completed PhCalculation
    """
    from aiida.orm.data.parameter import ParameterData
    from aiida.orm.utils import CalculationFactory, WorkflowFactory
    from aiida.work.run import run
    from aiida_quantumespresso.utils.resources import get_default_options

    Q2rBaseWorkChain = WorkflowFactory('quantumespresso.q2r.base')

    options = get_default_options(max_num_machines, max_wallclock_seconds)

    inputs = {
        'code': code,
        'parent_folder': parent_calc.out.retrieved,
        'options': ParameterData(dict=options),
    }

    run(Q2rBaseWorkChain, **inputs)