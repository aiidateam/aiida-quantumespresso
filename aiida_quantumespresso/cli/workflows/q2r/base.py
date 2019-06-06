# -*- coding: utf-8 -*-
"""Command line scripts to launch a `Q2rBaseWorkChain` for testing and demonstration purposes."""
from __future__ import absolute_import

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from ...cli import workflow_launch
from ...utils import launch
from ...utils import options as options_qe


@workflow_launch.command('q2r-base')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.q2r'))
@options.CALCULATION(type=types.CalculationParamType(sub_classes=('aiida.calculations:quantumespresso.ph',)))
@options_qe.CLEAN_WORKDIR()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def launch_workflow(code, calculation, clean_workdir, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run the `Q2rBaseWorkChain` for a previously completed `PhCalculation`."""
    from aiida.orm import Bool, Dict
    from aiida.plugins import WorkflowFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    inputs = {
        'code': code,
        'parent_folder': calculation.outputs.retrieved,
        'options': Dict(dict=get_default_options(max_num_machines, max_wallclock_seconds, with_mpi)),
    }

    if clean_workdir:
        inputs['clean_workdir'] = Bool(True)

    launch.launch_process(WorkflowFactory('quantumespresso.q2r.base'), daemon, **inputs)
