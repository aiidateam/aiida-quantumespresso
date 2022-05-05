# -*- coding: utf-8 -*-
"""Command line scripts to launch a `Q2rBaseWorkChain` for testing and demonstration purposes."""
from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators
import click

from .. import cmd_launch
from ...utils import launch
from ...utils import options as options_qe


@cmd_launch.command('q2r-base')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.q2r'))
@options.CALCULATION(required=True)
@options_qe.CLEAN_WORKDIR()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def launch_workflow(code, calculation, clean_workdir, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run the `Q2rBaseWorkChain` for a previously completed `PhCalculation`."""
    from aiida.orm import Bool
    from aiida.plugins import WorkflowFactory

    from aiida_quantumespresso.utils.resources import get_default_options

    expected_process_type = 'aiida.calculations:quantumespresso.ph'
    if calculation.process_type != expected_process_type:
        raise click.BadParameter(
            f'input calculation node has process_type: {calculation.process_type}; should be {expected_process_type}'
        )

    inputs = {
        'q2r': {
            'code': code,
            'parent_folder': calculation.outputs.remote_folder,
            'metadata': {
                'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
            }
        }
    }

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    launch.launch_process(WorkflowFactory('quantumespresso.q2r.base'), daemon, **inputs)
