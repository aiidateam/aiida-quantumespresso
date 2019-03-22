# -*- coding: utf-8 -*-
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
def launch(
    code, calculation, kpoints_mesh, clean_workdir, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """
    Run the MatdynBaseWorkChain for a previously completed Q2rCalculation
    """
    from aiida.orm.nodes.data.base import Bool
    from aiida.orm.nodes.data.dict import Dict
    from aiida.plugins import WorkflowFactory
    from aiida.work import launch
    from aiida_quantumespresso.utils.resources import get_default_options

    MatdynBaseWorkChain = WorkflowFactory('quantumespresso.matdyn.base')

    options = get_default_options(max_num_machines, max_wallclock_seconds)

    inputs = {
        'code': code,
        'kpoints': kpoints_mesh,
        'parent_folder': calculation.out.force_constants,
        'options': Dict(dict=options),
    }

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    if daemon:
        workchain = launch.submit(MatdynBaseWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(MatdynBaseWorkChain.__name__, workchain.pk))
    else:
        launch.run(MatdynBaseWorkChain, **inputs)
