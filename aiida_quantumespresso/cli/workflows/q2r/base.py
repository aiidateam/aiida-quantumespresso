# -*- coding: utf-8 -*-
"""Command line scripts to launch a `Q2rBaseWorkChain` for testing and demonstration purposes."""
from __future__ import absolute_import
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
def cli(code, calculation, clean_workdir, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run the `Q2rBaseWorkChain` for a previously completed `PhCalculation`."""
    from aiida.engine import launch
    from aiida.orm import Bool, Dict
    from aiida.plugins import WorkflowFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    Q2rBaseWorkChain = WorkflowFactory('quantumespresso.q2r.base')  # pylint: disable=invalid-name

    inputs = {
        'code': code,
        'parent_folder': calculation.out.retrieved,
        'options': Dict(dict=get_default_options(max_num_machines, max_wallclock_seconds, with_mpi)),
    }

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    if daemon:
        workchain = launch.submit(Q2rBaseWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(Q2rBaseWorkChain.__name__, workchain.pk))
    else:
        launch.run(Q2rBaseWorkChain, **inputs)
