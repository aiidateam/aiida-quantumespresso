# -*- coding: utf-8 -*-
"""Command line scripts to launch a `MatdynBaseWorkChain` for testing and demonstration purposes."""
from __future__ import absolute_import
import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from aiida_quantumespresso.cli.utils import options as options_qe


@click.command()
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.matdyn'))
@options.CALCULATION(type=types.CalculationParamType(sub_classes=('aiida.calculations:quantumespresso.q2r',)))
@options_qe.KPOINTS_MESH(default=[2, 2, 2])
@options_qe.CLEAN_WORKDIR()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def cli(code, calculation, kpoints_mesh, clean_workdir, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run the `MatdynBaseWorkChain` for a previously completed `Q2rCalculation`."""
    from aiida.engine import launch
    from aiida.orm import Bool, Dict
    from aiida.plugins import WorkflowFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    MatdynBaseWorkChain = WorkflowFactory('quantumespresso.matdyn.base')  # pylint: disable=invalid-name

    inputs = {
        'code': code,
        'kpoints': kpoints_mesh,
        'parent_folder': calculation.outputs.force_constants,
        'options': Dict(dict=get_default_options(max_num_machines, max_wallclock_seconds, with_mpi)),
    }

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    if daemon:
        workchain = launch.submit(MatdynBaseWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(MatdynBaseWorkChain.__name__, workchain.pk))
    else:
        launch.run(MatdynBaseWorkChain, **inputs)
