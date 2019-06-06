# -*- coding: utf-8 -*-
"""Command line scripts to launch a `MatdynCalculation` for testing and demonstration purposes."""
from __future__ import absolute_import

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from ..cli import calculation_launch
from ..utils import launch
from ..utils import options as options_qe


@calculation_launch.command('matdyn')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.matdyn'))
@options.DATUM(
    required=True,
    type=types.DataParamType(sub_classes=('aiida.data:quantumespresso.forceconstants',)),
    help='A ForceconstantsData node produced by a `Q2rCalculation`')
@options_qe.KPOINTS_MESH(default=[1, 1, 1])
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def launch_calculation(code, datum, kpoints_mesh, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run a MatdynCalculation."""
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    inputs = {
        'code': code,
        'kpoints': kpoints_mesh,
        'parent_folder': datum,
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    launch.launch_process(CalculationFactory('quantumespresso.matdyn'), daemon, **inputs)
