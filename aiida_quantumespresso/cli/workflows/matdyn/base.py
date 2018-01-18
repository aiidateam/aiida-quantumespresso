# -*- coding: utf-8 -*-
import click
from aiida_quantumespresso.utils.click import command
from aiida_quantumespresso.utils.click import options


@command()
@options.code()
@options.parent_calc(callback_kwargs={'entry_point': 'quantumespresso.q2r'})
@options.kpoint_mesh()
@options.max_num_machines()
@options.max_wallclock_seconds()
def launch(
    code, parent_calc, kpoints, max_num_machines, max_wallclock_seconds):
    """
    Run the MatdynBaseWorkChain for a previously completed Q2rCalculation
    """
    from aiida.orm.data.parameter import ParameterData
    from aiida.orm.utils import CalculationFactory, WorkflowFactory
    from aiida.work.run import run
    from aiida_quantumespresso.utils.resources import get_default_options

    MatdynBaseWorkChain = WorkflowFactory('quantumespresso.matdyn.base')

    options = get_default_options(max_num_machines, max_wallclock_seconds)

    inputs = {
        'code': code,
        'kpoints': kpoints,
        'parent_folder': parent_calc.out.force_constants,
        'options': ParameterData(dict=options),
    }

    run(MatdynBaseWorkChain, **inputs)