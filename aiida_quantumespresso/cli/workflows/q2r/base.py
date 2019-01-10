# -*- coding: utf-8 -*-
import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from aiida_quantumespresso.cli.utils import options as options_qe


@click.command()
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.q2r'))
@options.CALCULATION(type=types.CalculationParamType(sub_classes=('aiida.calculations:quantumespresso.ph',)))
@options_qe.CLEAN_WORKDIR()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def launch(
    code, calculation, clean_workdir, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """
    Run the Q2rBaseWorkChain for a previously completed PhCalculation
    """
    from aiida.orm.data.base import Bool
    from aiida.orm.data.parameter import ParameterData
    from aiida.orm.utils import WorkflowFactory
    from aiida.work import launch
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
        workchain = launch.submit(Q2rBaseWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(Q2rBaseWorkChain.__name__, workchain.pk))
    else:
        launch.run(Q2rBaseWorkChain, **inputs)
