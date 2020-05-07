# -*- coding: utf-8 -*-
"""Command line scripts to launch a `NebCalculation` for testing and demonstration purposes."""
from __future__ import absolute_import

import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from ..utils import launch
from ..utils import options as options_qe
from ..utils import validate
from . import cmd_launch


@cmd_launch.command('neb')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.neb'))
@click.option(
    '-s',
    '--structures',
    nargs=2,
    type=types.DataParamType(sub_classes=('aiida.data:structure',)),
    help='Two StructureData nodes representing the initial and final structures',
    metavar='<FIRST LAST>',
    required=True,
)
@click.option(
    '-I',
    '--num-images',
    type=click.INT,
    help='Number of points (images) used to discretize the path',
    default=3,
)
@click.option(
    '-N',
    '--num-steps',
    type=click.INT,
    help='Maximum number of path optimization steps',
    default=20,
)
@options_qe.PSEUDO_FAMILY(required=True)
@options_qe.KPOINTS_MESH(default=[2, 2, 2])
@options_qe.ECUTWFC()
@options_qe.ECUTRHO()
@options_qe.SMEARING()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@options_qe.PARENT_FOLDER()
@options.DRY_RUN()
@decorators.with_dbenv()
def launch_calculation(
    code, structures, num_images, num_steps, pseudo_family, kpoints_mesh, ecutwfc, ecutrho, smearing, max_num_machines,
    max_wallclock_seconds, with_mpi, daemon, parent_folder, dry_run
):
    """Run a NebCalculation.

    Note that some parameters are hardcoded.
    """
    from aiida.orm import Dict
    from aiida.orm.nodes.data.upf import get_pseudos_from_structure
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    pw_parameters = {
        'CONTROL': {
            'calculation': 'relax',
        },
        'SYSTEM': {
            'ecutwfc': ecutwfc,
            'ecutrho': ecutrho,
        }
    }

    neb_parameters = {
        'PATH': {
            'restart_mode': ('restart' if parent_folder else 'from_scratch'),
            'opt_scheme': 'broyden',
            'num_of_images': num_images,
            'nstep_path': num_steps,
        },
    }

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
            # would be confusing ("Submitted NebCalculation<None> to the daemon")
            raise click.BadParameter('cannot send to the daemon if in dry_run mode', param_hint='--daemon')
        inputs['metadata']['store_provenance'] = False
        inputs['metadata']['dry_run'] = True

    launch.launch_process(CalculationFactory('quantumespresso.neb'), daemon, **inputs)
