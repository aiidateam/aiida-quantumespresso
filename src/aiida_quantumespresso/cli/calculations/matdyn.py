# -*- coding: utf-8 -*-
"""Command line scripts to launch a `MatdynCalculation` for testing and demonstration purposes."""
from aiida.cmdline.params import options as options_core
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators

from . import cmd_launch
from ..utils import launch, options


@cmd_launch.command('matdyn')
@options_core.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.matdyn'))
@options_core.DATUM(
    required=True,
    type=types.DataParamType(sub_classes=('aiida.data:quantumespresso.force_constants',)),
    help='A ForceConstantsData node produced by a `Q2rCalculation`'
)
@options.KPOINTS_MESH(default=[1, 1, 1])
@options.MAX_NUM_MACHINES()
@options.MAX_WALLCLOCK_SECONDS()
@options.WITH_MPI()
@options.DAEMON()
@decorators.with_dbenv()
def launch_calculation(code, datum, kpoints_mesh, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run a MatdynCalculation."""
    from aiida.plugins import CalculationFactory

    from aiida_quantumespresso.utils.resources import get_default_options

    inputs = {
        'code': code,
        'kpoints': kpoints_mesh,
        'force_constants': datum,
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    launch.launch_process(CalculationFactory('quantumespresso.matdyn'), daemon, **inputs)
