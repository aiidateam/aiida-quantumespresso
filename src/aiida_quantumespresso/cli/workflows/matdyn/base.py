# -*- coding: utf-8 -*-
"""Command line scripts to launch a `MatdynBaseWorkChain` for testing and demonstration purposes."""
from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from .. import cmd_launch
from ...utils import launch
from ...utils import options as options_qe


@cmd_launch.command('matdyn-base')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.matdyn'))
@options.DATUM(
    required=True,
    type=types.DataParamType(sub_classes=('aiida.data:quantumespresso.force_constants',)),
    help='A ForceConstantsData node produced by a `Q2rCalculation`'
)
@options_qe.KPOINTS_MESH(default=[2, 2, 2])
@options_qe.CLEAN_WORKDIR()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def launch_workflow(
    code, datum, kpoints_mesh, clean_workdir, max_num_machines, max_wallclock_seconds, with_mpi, daemon
):
    """Run the `MatdynBaseWorkChain` for a previously completed `Q2rCalculation`."""
    from aiida.orm import Bool
    from aiida.plugins import WorkflowFactory

    from aiida_quantumespresso.utils.resources import get_default_options

    inputs = {
        'matdyn': {
            'code': code,
            'kpoints': kpoints_mesh,
            'force_constants': datum,
            'metadata': {
                'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
            }
        }
    }

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    launch.launch_process(WorkflowFactory('quantumespresso.matdyn.base'), daemon, **inputs)
