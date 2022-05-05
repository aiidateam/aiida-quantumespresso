# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PhBaseWorkChain` for testing and demonstration purposes."""
from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from .. import cmd_launch
from ...utils import launch
from ...utils import options as options_qe


@cmd_launch.command('ph-base')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.ph'))
@options.CALCULATION(required=True)
@options_qe.KPOINTS_MESH(default=[2, 2, 2])
@options_qe.CLEAN_WORKDIR()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def launch_workflow(
    code, calculation, kpoints_mesh, clean_workdir, max_num_machines, max_wallclock_seconds, with_mpi, daemon
):
    """Run the `PhBaseWorkChain` for a previously completed `PwCalculation`."""
    from aiida.orm import Bool, Dict
    from aiida.plugins import WorkflowFactory

    from aiida_quantumespresso.utils.resources import get_default_options

    inputs = {
        'ph': {
            'code': code,
            'qpoints': kpoints_mesh,
            'parent_folder': calculation.outputs.remote_folder,
            'parameters': Dict({'INPUTPH': {}}),
            'metadata': {
                'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
            }
        }
    }

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    launch.launch_process(WorkflowFactory('quantumespresso.ph.base'), daemon, **inputs)
