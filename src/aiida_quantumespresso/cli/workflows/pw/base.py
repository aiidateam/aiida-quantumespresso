# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PwBaseWorkChain` for testing and demonstration purposes."""
from aiida.cmdline.params import options as options_core
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators
import click

from .. import cmd_launch
from ...utils import defaults, launch, options, validate


@cmd_launch.command('pw-base')
@options_core.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.pw'))
@options.STRUCTURE(default=defaults.get_structure)
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
@decorators.with_dbenv()
def launch_workflow(
    code, structure, pseudo_family, kpoints_distance, ecutwfc, ecutrho, hubbard_u, hubbard_v, hubbard_file_pk,
    starting_magnetization, smearing, automatic_parallelization, clean_workdir, max_num_machines, max_wallclock_seconds,
    with_mpi, daemon
):
    """Run a `PwBaseWorkChain`."""
    from aiida.orm import Bool, Dict, Float
    from aiida.plugins import WorkflowFactory

    from aiida_quantumespresso.utils.resources import get_automatic_parallelization_options, get_default_options

    builder = WorkflowFactory('quantumespresso.pw.base').get_builder()

    cutoff_wfc, cutoff_rho = pseudo_family.get_recommended_cutoffs(structure=structure, unit='Ry')

    parameters = {
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

    builder.pw.code = code
    builder.pw.structure = structure
    builder.pw.parameters = Dict(parameters)
    builder.pw.pseudos = pseudo_family.get_pseudos(structure=structure)
    builder.kpoints_distance = Float(kpoints_distance)

    if hubbard_file:
        builder.hubbard_file = hubbard_file

    if automatic_parallelization:
        automatic_parallelization = get_automatic_parallelization_options(max_num_machines, max_wallclock_seconds)
        builder.automatic_parallelization = Dict(automatic_parallelization)
    else:
        builder.pw.metadata.options = get_default_options(max_num_machines, max_wallclock_seconds, with_mpi)

    if clean_workdir:
        builder.clean_workdir = Bool(True)

    launch.launch_process(builder, daemon)
