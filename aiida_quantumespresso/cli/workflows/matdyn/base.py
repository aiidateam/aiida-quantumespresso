# -*- coding: utf-8 -*-
import click
from aiida.utils.cli import command
from aiida.utils.cli import options
from aiida_quantumespresso.utils.cli import options as options_qe


@command()
@options.code(callback_kwargs={'entry_point': 'quantumespresso.matdyn'})
@options.calculation(callback_kwargs={'entry_point': 'quantumespresso.q2r'})
@options.kpoint_mesh()
@options.max_num_machines()
@options.max_wallclock_seconds()
@options.daemon()
@options_qe.clean_workdir()
def launch(
    code, calculation, kpoints, max_num_machines, max_wallclock_seconds, daemon, clean_workdir):
    """
    Run the MatdynBaseWorkChain for a previously completed Q2rCalculation
    """
    from aiida.orm.data.base import Bool
    from aiida.orm.data.parameter import ParameterData
    from aiida.orm.utils import CalculationFactory, WorkflowFactory
    from aiida.work.launch import run, submit
    from aiida_quantumespresso.utils.resources import get_default_options

    MatdynBaseWorkChain = WorkflowFactory('quantumespresso.matdyn.base')

    options = get_default_options(max_num_machines, max_wallclock_seconds)

    inputs = {
        'code': code,
        'kpoints': kpoints,
        'parent_folder': calculation.out.force_constants,
        'options': ParameterData(dict=options),
    }

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    if daemon:
        workchain = submit(MatdynBaseWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(MatdynBaseWorkChain.__name__, workchain.pk))
    else:
        run(MatdynBaseWorkChain, **inputs)
