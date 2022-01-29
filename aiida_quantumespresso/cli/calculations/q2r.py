# -*- coding: utf-8 -*-
"""Command line scripts to launch a `Q2rCalculation` for testing and demonstration purposes."""
from aiida.cmdline.params import options as options_core
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators
import click

from . import cmd_launch
from ..utils import launch, options


@cmd_launch.command('q2r')
@options_core.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.q2r'))
@options_core.CALCULATION(required=True)
@options.MAX_NUM_MACHINES()
@options.MAX_WALLCLOCK_SECONDS()
@options.WITH_MPI()
@options.DAEMON()
@decorators.with_dbenv()
def launch_calculation(code, calculation, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run a Q2rCalculation."""
    from aiida.plugins import CalculationFactory

    from aiida_quantumespresso.utils.resources import get_default_options

    # Check that the parent calculation node comes from quantumespresso.ph.
    # I cannot move this check into the option declaration, because CalcJobNode is not subclassed by the specific
    # calculation plugins (only Process is), and there is no feature yet to filter by the associated process_type.
    expected_process_type = 'aiida.calculations:quantumespresso.ph'
    if calculation.process_type != expected_process_type:
        raise click.BadParameter(
            f'input calculation node has process_type: {calculation.process_type}; should be {expected_process_type}'
        )

    inputs = {
        'code': code,
        'parent_folder': calculation.outputs.remote_folder,
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    launch.launch_process(CalculationFactory('quantumespresso.q2r'), daemon, **inputs)
