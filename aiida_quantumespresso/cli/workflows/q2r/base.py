# -*- coding: utf-8 -*-
import click
from aiida.utils.cli import command
from aiida.utils.cli import options
from aiida_quantumespresso.utils.cli import options as options_qe



@command()
@options.code(callback_kwargs={'entry_point': 'quantumespresso.q2r'})
@options.calculation(callback_kwargs={'entry_point': 'quantumespresso.ph'})
@options.max_num_machines()
@options.max_wallclock_seconds()
@options.daemon()
@options_qe.clean_workdir()
def launch(
    code, calculation, max_num_machines, max_wallclock_seconds, daemon, clean_workdir):
    """
    Run the Q2rBaseWorkChain for a previously completed PhCalculation
    """
    from aiida.orm.data.base import Bool
    from aiida.orm.data.parameter import ParameterData
    from aiida.orm.utils import CalculationFactory, WorkflowFactory
    from aiida.work.launch import run, submit
    from aiida_quantumespresso.utils.resources import get_default_options

    Q2rBaseWorkChain = WorkflowFactory('quantumespresso.q2r.base')

    options = get_default_options(max_num_machines, max_wallclock_seconds)

    inputs = {
        'code': code,
        'parent_folder': calculation.out.retrieved,
        'options': ParameterData(dict=options),
    }

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    if daemon:
        workchain = submit(Q2rBaseWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(Q2rBaseWorkChain.__name__, workchain.pk))
    else:
        run(Q2rBaseWorkChain, **inputs)
