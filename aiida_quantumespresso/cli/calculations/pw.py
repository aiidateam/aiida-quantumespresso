# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PwCalculation` for testing and demonstration purposes."""
from __future__ import absolute_import

import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from ..utils import defaults
from ..utils import launch
from ..utils import options as options_qe
from ..utils import validate
from . import cmd_launch

CALCS_REQUIRING_PARENT = set(['nscf'])


@cmd_launch.command('pw')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.pw'))
@options_qe.STRUCTURE(default=defaults.get_structure)
@options_qe.PSEUDO_FAMILY(required=True)
@options_qe.KPOINTS_MESH(default=[2, 2, 2])
@options_qe.ECUTWFC()
@options_qe.ECUTRHO()
@options_qe.HUBBARD_U()
@options_qe.HUBBARD_V()
@options_qe.HUBBARD_FILE()
@options_qe.STARTING_MAGNETIZATION()
@options_qe.SMEARING()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@options_qe.PARENT_FOLDER()
@options.DRY_RUN()
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
    from aiida.orm import Dict
    from aiida.orm.nodes.data.upf import get_pseudos_from_structure
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    parameters = {
        'CONTROL': {
            'calculation': mode,
        },
        'SYSTEM': {
            'ecutwfc': ecutwfc,
            'ecutrho': ecutrho,
        }
    }

    if mode in CALCS_REQUIRING_PARENT and not parent_folder:
        raise click.BadParameter("calculation '{}' requires a parent folder".format(mode), param_hint='--parent-folder')

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
        from aiida.orm import KpointsData
        unfolded_list = kpoints_mesh.get_kpoints_mesh(print_list=True)
        kpoints_mesh = KpointsData()
        kpoints_mesh.set_kpoints(unfolded_list)

    inputs = {
        'code': code,
        'structure': structure,
        'pseudos': get_pseudos_from_structure(structure, pseudo_family),
        'kpoints': kpoints_mesh,
        'parameters': Dict(dict=parameters),
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
