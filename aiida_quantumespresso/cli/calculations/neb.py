# -*- coding: utf-8 -*-
"""Command line scripts to launch a `NebCalculation` for testing and demonstration purposes."""
from __future__ import absolute_import

import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from ..cli import calculation_launch
from ..utils import launch
from ..utils import options as options_qe
from ..utils import validate


@calculation_launch.command('neb')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.neb'))
@click.option(
    '-s',
    '--structures',
    nargs=2,
    type=types.DataParamType(sub_classes=('aiida.data:structure',)),
    help='Two StructureData nodes, the initial and final structures',
    required=True
)
@options_qe.PSEUDO_FAMILY(required=True)
@options_qe.KPOINTS_MESH(default=[2, 2, 2])
@options_qe.ECUTWFC()
@options_qe.ECUTRHO()
@options_qe.STARTING_MAGNETIZATION()
@options_qe.SMEARING()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@options_qe.PARENT_FOLDER()
@options.DRY_RUN()
@decorators.with_dbenv()
def launch_calculation(
    code, structures, pseudo_family, kpoints_mesh, ecutwfc, ecutrho, starting_magnetization, smearing, max_num_machines,
    max_wallclock_seconds, with_mpi, daemon, parent_folder, dry_run
):
    """
    Run a NebCalculation.
    Note that some parameters are hardcoded.
    """
    from aiida.orm import Dict
    from aiida.orm.nodes.data.upf import get_pseudos_from_structure
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    # TODO: move some hardcoded parameters to unit tests
    # TODO: because I cannot give a structure-specific settings dict, I can't apply these constraints only to one of the
    #  two boundary images. Do I need to?
    settings = {
        'fixed_coords':
            [[False, True, True],
             [True, True, True],
             [False, True, True]],
        'CLIMBING_IMAGES': [4],
    }

    pw_parameters = {
        'CONTROL': {
            'calculation': 'relax',  # TODO: this needs to be set automatically by the calculation class (is relax correct?)
        },
        'SYSTEM': {
            'ecutwfc': ecutwfc,
            'ecutrho': ecutrho,
            'occupations': 'smearing',
            'degauss': 0.003,
            'nspin': 2,
            'starting_magnetization': 0.5,
        },
        'ELECTRONS': {
            'conv_thr': 1e-8,
            'mixing_beta': 0.3,
        }
    }

    neb_parameters = {
        'PATH': {
            'restart_mode': ('restart' if parent_folder else 'from_scratch'),
            'nstep_path': 20,
            'ds': 2.,
            'opt_scheme': 'broyden',
            'num_of_images': 6,
            'k_max': 0.3,
            'k_min': 0.2,
            'CI_scheme': 'manual',  # 'manual', 'auto', or 'no-CI'
            'path_thr': 0.05,
        },
    }

    for structure in structures:
        try:
            validate.validate_starting_magnetization(structure, pw_parameters, starting_magnetization)
        except ValueError as exception:
            raise click.BadParameter(str(exception))

        try:
            validate.validate_smearing(pw_parameters, smearing)
        except ValueError as exception:
            raise click.BadParameter(str(exception))

    inputs = {
        'code': code,
        'first_structure': structures[0],
        'last_structure': structures[1],
        'pw': {
            'pseudos': get_pseudos_from_structure(structures[0], pseudo_family),
            'kpoints': kpoints_mesh,
            'parameters': Dict(dict=pw_parameters),
        },
        'parameters': Dict(dict=neb_parameters),
        'settings': Dict(dict=settings),
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    if parent_folder:
        inputs['parent_folder'] = parent_folder

    if dry_run:
        if daemon:
            # .submit() would forward to .run(), but it's better to stop here,
            # since it's a bit unexpected and the log messages output to screen
            # would be confusing ("Submitted PwCalculation<None> to the daemon")
            raise click.BadParameter('cannot send to the daemon if in dry_run mode', param_hint='--daemon')
        inputs.setdefault('metadata', {})['store_provenance'] = False
        inputs['metadata']['dry_run'] = True

    launch.launch_process(CalculationFactory('quantumespresso.neb'), daemon, **inputs)
