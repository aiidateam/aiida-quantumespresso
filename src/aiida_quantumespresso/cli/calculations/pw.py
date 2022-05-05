# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PwCalculation` for testing and demonstration purposes."""
from aiida.cmdline.params import options as options_core
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators
import click

from . import cmd_launch
from ..utils import defaults, launch, options, validate

CALCS_REQUIRING_PARENT = set(['nscf'])


@cmd_launch.command('pw')
@options_core.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.pw'))
@options.STRUCTURE(default=defaults.get_structure)
@options.PSEUDO_FAMILY()
@options.KPOINTS_MESH(default=[2, 2, 2])
@options.ECUTWFC()
@options.ECUTRHO()
@options.HUBBARD_U()
@options.HUBBARD_V()
@options.HUBBARD_FILE()
@options.STARTING_MAGNETIZATION()
@options.SMEARING()
@options.MAX_NUM_MACHINES()
@options.MAX_WALLCLOCK_SECONDS()
@options.WITH_MPI()
@options.DAEMON()
@options.PARENT_FOLDER()
@options_core.DRY_RUN()
@click.option(
    '-z',
    '--calculation-mode',
    'mode',
    type=click.Choice(['scf', 'nscf', 'relax', 'vc-relax']),
    default='scf',
    show_default=True,
    help='Select the calculation mode.'
)
@click.option(
    '-u',
    '--unfolded-kpoints',
    'unfolded_kpoints',
    is_flag=True,
    help='Unfold the k-points grid to the whole grid without reducing it by symmetry (useful mainly for NSCF).'
)
@decorators.with_dbenv()
def launch_calculation(
    code, structure, pseudo_family, kpoints_mesh, ecutwfc, ecutrho, hubbard_u, hubbard_v, hubbard_file_pk,
    starting_magnetization, smearing, max_num_machines, max_wallclock_seconds, with_mpi, daemon, parent_folder, dry_run,
    mode, unfolded_kpoints
):
    """Run a PwCalculation."""
    from aiida.orm import Dict, KpointsData
    from aiida.plugins import CalculationFactory

    from aiida_quantumespresso.utils.resources import get_default_options

    cutoff_wfc, cutoff_rho = pseudo_family.get_recommended_cutoffs(structure=structure, unit='Ry')

    parameters = {
        'CONTROL': {
            'calculation': mode,
        },
        'SYSTEM': {
            'ecutwfc': ecutwfc or cutoff_wfc,
            'ecutrho': ecutrho or cutoff_rho,
        }
    }

    if mode in CALCS_REQUIRING_PARENT and not parent_folder:
        raise click.BadParameter(f"calculation '{mode}' requires a parent folder", param_hint='--parent-folder')

    try:
        hubbard_file = validate.validate_hubbard_parameters(
            structure, parameters, hubbard_u, hubbard_v, hubbard_file_pk
        )
    except ValueError as exception:
        raise click.BadParameter(str(exception))

    try:
        validate.validate_starting_magnetization(structure, parameters, starting_magnetization)
    except ValueError as exception:
        raise click.BadParameter(str(exception))

    try:
        validate.validate_smearing(parameters, smearing)
    except ValueError as exception:
        raise click.BadParameter(str(exception))

    if unfolded_kpoints:
        unfolded_list = kpoints_mesh.get_kpoints_mesh(print_list=True)
        kpoints_mesh = KpointsData()
        kpoints_mesh.set_kpoints(unfolded_list)

    inputs = {
        'code': code,
        'structure': structure,
        'pseudos': pseudo_family.get_pseudos(structure=structure),
        'kpoints': kpoints_mesh,
        'parameters': Dict(parameters),
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    if parent_folder:
        inputs['parent_folder'] = parent_folder

    if hubbard_file:
        inputs['hubbard_file'] = hubbard_file

    if dry_run:
        if daemon:
            # .submit() would forward to .run(), but it's better to stop here,
            # since it's a bit unexpected and the log messages output to screen
            # would be confusing ("Submitted PwCalculation<None> to the daemon")
            raise click.BadParameter('cannot send to the daemon if in dry_run mode', param_hint='--daemon')
        inputs['metadata']['store_provenance'] = False
        inputs['metadata']['dry_run'] = True

    launch.launch_process(CalculationFactory('quantumespresso.pw'), daemon, **inputs)
