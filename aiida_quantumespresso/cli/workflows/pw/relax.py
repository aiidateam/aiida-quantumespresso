# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PwRelaxWorkChain` for testing and demonstration purposes."""
import click

from aiida.cmdline.params import options as options_core
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators

from ...utils import launch
from ...utils import options
from ...utils import validate
from .. import cmd_launch


@cmd_launch.command('pw-relax')
@options_core.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.pw'))
@options.STRUCTURE(required=True)
@options.PSEUDO_FAMILY()
@options.KPOINTS_DISTANCE()
@options.ECUTWFC()
@options.ECUTRHO()
@options.HUBBARD_U()
@options.HUBBARD_V()
@options.HUBBARD_FILE()
@options.STARTING_MAGNETIZATION()
@options.SMEARING()
@options.AUTOMATIC_PARALLELIZATION()
@options.CLEAN_WORKDIR()
@options.MAX_NUM_MACHINES()
@options.MAX_WALLCLOCK_SECONDS()
@options.WITH_MPI()
@options.DAEMON()
@click.option(
    '-f',
    '--final-scf',
    is_flag=True,
    default=False,
    show_default=True,
    help='Run a final scf calculation for the final relaxed structure.'
)
@decorators.with_dbenv()
def launch_workflow(
    code, structure, pseudo_family, kpoints_distance, ecutwfc, ecutrho, hubbard_u, hubbard_v, hubbard_file_pk,
    starting_magnetization, smearing, automatic_parallelization, clean_workdir, max_num_machines, max_wallclock_seconds,
    with_mpi, daemon, final_scf
):
    """Run a `PwRelaxWorkChain`."""
    from aiida.orm import Bool, Float, Dict, Str
    from aiida.plugins import WorkflowFactory

    from aiida_quantumespresso.utils.resources import get_default_options, get_automatic_parallelization_options

    builder = WorkflowFactory('quantumespresso.pw.relax').get_builder()

    cutoff_wfc, cutoff_rho = pseudo_family.get_recommended_cutoffs(structure=structure, unit='Ry')

    parameters = {
        'CONTROL': {
            'calculation': 'relax',
        },
        'SYSTEM': {
            'ecutwfc': ecutwfc or cutoff_wfc,
            'ecutrho': ecutrho or cutoff_rho,
        },
    }

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

    builder.structure = structure
    builder.base.kpoints_distance = Float(kpoints_distance)
    builder.base.pw.code = code
    builder.base.pw.pseudos = pseudo_family.get_pseudos(structure=structure)
    builder.base.pw.parameters = Dict(dict=parameters)

    if hubbard_file:
        builder.base.pw.hubbard_file = hubbard_file

    if automatic_parallelization:
        automatic_parallelization = get_automatic_parallelization_options(max_num_machines, max_wallclock_seconds)
        builder.base.automatic_parallelization = Dict(dict=automatic_parallelization)
    else:
        builder.base.pw.metadata.options = get_default_options(max_num_machines, max_wallclock_seconds, with_mpi)

    if clean_workdir:
        builder.clean_workdir = Bool(True)

    if final_scf:
        builder.base_final_scf.pseudo_family = Str(pseudo_family)
        builder.base_final_scf.kpoints_distance = Float(kpoints_distance)
        builder.base_final_scf.pw.code = code
        builder.base_final_scf.pw.parameters = Dict(dict=parameters)
        builder.base_final_scf.pw.metadata.options = get_default_options(
            max_num_machines, max_wallclock_seconds, with_mpi
        )

    launch.launch_process(builder, daemon)
