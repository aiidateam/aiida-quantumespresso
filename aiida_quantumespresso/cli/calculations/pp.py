# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PpCalculation` for testing and demonstration purposes."""
import click

from aiida.cmdline.params import options as options_core
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators

from ..utils import launch
from ..utils import options
from . import cmd_launch


@cmd_launch.command('pp')
@options_core.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.pp'))
@options_core.CALCULATION(required=True)
@options.MAX_NUM_MACHINES()
@options.MAX_WALLCLOCK_SECONDS()
@options.WITH_MPI()
@options.DAEMON()
@decorators.with_dbenv()
def launch_calculation(code, calculation, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run a PpCalculation."""
    from aiida import orm
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    # Check that the parent calculation node comes from quantumespresso.pw.
    # I cannot move this check into the option declaration, because CalcJobNode is not subclassed by the specific
    # calculation plugins (only Process is), and there is no feature yet to filter by the associated process_type.
    expected_process_type = 'aiida.calculations:quantumespresso.pw'
    if calculation.process_type != expected_process_type:
        raise click.BadParameter(
            'The input calculation node has a process_type: {}; should be {}'.format(
                calculation.process_type, expected_process_type
            )
        )

    parent_folder = calculation.get_outgoing(node_class=orm.RemoteData, link_label_filter='remote_folder').one().node

    inputs = {
        'code': code,
        'parent_folder': parent_folder,
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    launch.launch_process(CalculationFactory('quantumespresso.pp'), daemon, **inputs)
