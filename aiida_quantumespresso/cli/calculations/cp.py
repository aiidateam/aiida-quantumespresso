# -*- coding: utf-8 -*-
"""Command line scripts to launch a `CpCalculation` for testing and demonstration purposes."""
from __future__ import absolute_import

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from ..cli import calculation_launch
from ..utils import launch
from ..utils import options as options_qe


def silicon_structure():
    """Return a `StructureData` representing bulk silicon.

    The database will first be queried for the existence of the structure in case it might have already been created
    before. If nothing is found, a new structure will be created

    :return: a `StructureData` representing bulk silicon
    """
    from ase.spacegroup import crystal
    from aiida.orm import QueryBuilder, StructureData

    label = '__example__{name}'.format(name=__name__)
    structure = QueryBuilder().append(StructureData, filters={'label': label}).first()
    if not structure:
        alat = 5.4
        ase_structure = crystal(
            'Si',
            [(0, 0, 0)],
            spacegroup=227,
            cellpar=[alat, alat, alat, 90, 90, 90],
            primitive_cell=True,
        )
        structure = StructureData(ase=ase_structure)
        structure.label = label
        structure.store()
    else:
        # .first() returns none if it doesn't exist or a list of one item if objects exist
        structure = structure[0]

    return structure.uuid


@calculation_launch.command('cp')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.cp'))
@options_qe.STRUCTURE(default=silicon_structure)
@options_qe.PSEUDO_FAMILY(required=True)
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def launch_calculation(code, structure, pseudo_family, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run a CpCalculation."""
    from aiida.orm import Dict
    from aiida.orm.nodes.data.upf import get_pseudos_from_structure
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    parameters = {
        'CONTROL': {
            'calculation': 'cp',
            'restart_mode': 'from_scratch',
            'wf_collect': False,
            'iprint': 1,
            'isave': 100,
            'dt': 3.0,
            'max_seconds': 25 * 60,
            'nstep': 10,
        },
        'SYSTEM': {
            'ecutwfc': 30.0,
            'ecutrho': 240.0,
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
        'pseudos': get_pseudos_from_structure(structure, pseudo_family),
        'parameters': Dict(dict=parameters),
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    launch.launch_process(CalculationFactory('quantumespresso.cp'), daemon, **inputs)
