# -*- coding: utf-8 -*-
import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from aiida_quantumespresso.cli.utils import options as options_qe


@click.command()
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.ph'))
@options.CALCULATION(type=types.CalculationParamType(sub_classes=('aiida.calculations:quantumespresso.pw',)))
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
    Run the PhBaseWorkChain for a previously completed PwCalculation
    """
    from aiida.orm.data.base import Bool
    from aiida.orm.data.parameter import ParameterData
    from aiida.orm.utils import WorkflowFactory
    from aiida.work import launch
    from aiida_quantumespresso.utils.resources import get_default_options

    PhBaseWorkChain = WorkflowFactory('quantumespresso.ph.base')

    parameters = {
        'INPUTPH': {
        }
    }

    options = get_default_options(max_num_machines, max_wallclock_seconds, with_mpi)

    inputs = {
        'code': code,
        'qpoints': kpoints_mesh,
        'parent_folder': calculation.out.remote_folder,
        'parameters': ParameterData(dict=parameters),
        'options': ParameterData(dict=options),
    }

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    if daemon:
        workchain = launch.submit(PhBaseWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(PhBaseWorkChain.__name__, workchain.pk))
    else:
        launch.run(PhBaseWorkChain, **inputs)
