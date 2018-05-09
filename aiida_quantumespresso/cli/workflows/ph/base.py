# -*- coding: utf-8 -*-
import click
from aiida_quantumespresso.utils.click import command
from aiida_quantumespresso.utils.click import options


@command()
@options.code()
@options.parent_calc(callback_kwargs={'entry_point': 'quantumespresso.pw'})
@options.kpoint_mesh()
@options.max_num_machines()
@options.max_wallclock_seconds()
@options.daemon()
def launch(
    code, parent_calc, kpoints, max_num_machines, max_wallclock_seconds, daemon):
    """
    Run the PhBaseWorkChain for a previously completed PwCalculation
    """
    from aiida.orm.data.parameter import ParameterData
    from aiida.orm.utils import CalculationFactory, WorkflowFactory
    from aiida.work.run import run, submit
    from aiida_quantumespresso.utils.resources import get_default_options

    PwCalculation = CalculationFactory('quantumespresso.pw')
    PhBaseWorkChain = WorkflowFactory('quantumespresso.ph.base')

    parameters = {
        'INPUTPH': {
        }
    }

    options = get_default_options(max_num_machines, max_wallclock_seconds)

    inputs = {
        'code': code,
        'qpoints': kpoints,
        'parent_folder': parent_calc.out.remote_folder,
        'parameters': ParameterData(dict=parameters),
        'options': ParameterData(dict=options),
    }

    if daemon:
        workchain = submit(PhBaseWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(PhBaseWorkChain.__name__, workchain.pid))
    else:
        run(PhBaseWorkChain, **inputs)
