# -*- coding: utf-8 -*-
"""Command line scripts to launch a `CpCalculation` for testing and demonstration purposes."""
from aiida.cmdline.params import options as options_core
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators

from . import cmd_launch
from ..utils import defaults, launch, options


@cmd_launch.command('cp')
@options_core.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.cp'))
@options.STRUCTURE(default=defaults.get_structure)
@options.PSEUDO_FAMILY()
@options.MAX_NUM_MACHINES()
@options.MAX_WALLCLOCK_SECONDS()
@options.WITH_MPI()
@options.DAEMON()
@decorators.with_dbenv()
def launch_calculation(code, structure, pseudo_family, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run a CpCalculation."""
    from aiida.orm import Dict
    from aiida.plugins import CalculationFactory

    from aiida_quantumespresso.utils.resources import get_default_options

    cutoff_wfc, cutoff_rho = pseudo_family.get_recommended_cutoffs(structure=structure, unit='Ry')

    parameters = {
        'CONTROL': {
            'calculation': 'cp',
            'restart_mode': 'from_scratch',
            'wf_collect': False,
            'iprint': 1,
            'isave': 100,
            'dt': 3.0,
            'nstep': 10,
        },
        'SYSTEM': {
            'ecutwfc': cutoff_wfc,
            'ecutrho': cutoff_rho,
            'nr1b': 24,
            'nr2b': 24,
            'nr3b': 24,
        },
        'ELECTRONS': {
            'electron_damping': 1.0e-1,
            'electron_dynamics': 'damp',
            'emass': 400.0,
            'emass_cutoff': 3.0,
        },
        'IONS': {
            'ion_dynamics': 'none'
        },
    }

    inputs = {
        'code': code,
        'structure': structure,
        'pseudos': pseudo_family.get_pseudos(structure=structure),
        'parameters': Dict(parameters),
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    launch.launch_process(CalculationFactory('quantumespresso.cp'), daemon, **inputs)
