# -*- coding: utf-8 -*-
"""Command line scripts to launch a `NebCalculation` for testing and demonstration purposes."""
from aiida.cmdline.params import options as options_core
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators
import click

from . import cmd_launch
from ..utils import launch, options, validate


@cmd_launch.command('neb')
@options_core.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.neb'))
@click.option(
    '-s',
    '--structures',
    nargs=2,
    type=types.DataParamType(sub_classes=('aiida.data:core.structure',)),
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
@options.PSEUDO_FAMILY()
@options.KPOINTS_MESH(default=[2, 2, 2])
@options.ECUTWFC()
@options.ECUTRHO()
@options.SMEARING()
@options.MAX_NUM_MACHINES()
@options.MAX_WALLCLOCK_SECONDS()
@options.WITH_MPI()
@options.DAEMON()
@options.PARENT_FOLDER()
@options_core.DRY_RUN()
@decorators.with_dbenv()
def launch_calculation(
    code, structures, num_images, num_steps, pseudo_family, kpoints_mesh, ecutwfc, ecutrho, smearing, max_num_machines,
    max_wallclock_seconds, with_mpi, daemon, parent_folder, dry_run
):
    """Run a NebCalculation.

    Note that some parameters are hardcoded.
    """
    from aiida.orm import Dict
    from aiida.plugins import CalculationFactory

    from aiida_quantumespresso.utils.resources import get_default_options

    cutoff_wfc, cutoff_rho = pseudo_family.get_recommended_cutoffs(structure=structures[0], unit='Ry')

    pw_parameters = {
        'CONTROL': {
            'calculation': 'relax',
        },
        'SYSTEM': {
            'ecutwfc': ecutwfc or cutoff_wfc,
            'ecutrho': ecutrho or cutoff_rho
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
            'pseudos': pseudo_family.get_pseudos(structure=structures[0]),
            'kpoints': kpoints_mesh,
            'parameters': Dict(pw_parameters),
        },
        'parameters': Dict(neb_parameters),
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
