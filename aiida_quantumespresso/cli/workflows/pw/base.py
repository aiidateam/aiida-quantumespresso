# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PwBaseWorkChain` for testing and demonstration purposes."""
from __future__ import absolute_import
import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from ...utils import defaults
from ...utils import launch
from ...utils import options as options_qe
from ...utils import validate
from .. import cmd_launch


@cmd_launch.command('pw-base')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.pw'))
@options_qe.STRUCTURE(default=defaults.get_structure)
@options_qe.PSEUDO_FAMILY(required=True)
@options_qe.KPOINTS_DISTANCE()
@options_qe.ECUTWFC()
@options_qe.ECUTRHO()
@options_qe.HUBBARD_U()
@options_qe.HUBBARD_V()
@options_qe.HUBBARD_FILE()
@options_qe.STARTING_MAGNETIZATION()
@options_qe.SMEARING()
@options_qe.AUTOMATIC_PARALLELIZATION()
@options_qe.CLEAN_WORKDIR()
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def launch_workflow(
    code, structure, pseudo_family, kpoints_distance, ecutwfc, ecutrho, hubbard_u, hubbard_v, hubbard_file_pk,
    starting_magnetization, smearing, automatic_parallelization, clean_workdir, max_num_machines, max_wallclock_seconds,
    with_mpi, daemon
):
    """Run a `PwBaseWorkChain`."""
    from aiida.orm import Bool, Float, Str, Dict
    from aiida.plugins import WorkflowFactory
    from aiida_quantumespresso.utils.resources import get_default_options, get_automatic_parallelization_options

    builder = WorkflowFactory('quantumespresso.pw.base').get_builder()

    parameters = {
        'SYSTEM': {
            'ecutwfc': ecutwfc,
            'ecutrho': ecutrho,
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
    builder.pw.parameters = Dict(dict=parameters)
    builder.pseudo_family = Str(pseudo_family)
    builder.kpoints_distance = Float(kpoints_distance)

    if hubbard_file:
        builder.hubbard_file = hubbard_file

    if automatic_parallelization:
        automatic_parallelization = get_automatic_parallelization_options(max_num_machines, max_wallclock_seconds)
        builder.automatic_parallelization = Dict(dict=automatic_parallelization)
    else:
        builder.pw.metadata.options = get_default_options(max_num_machines, max_wallclock_seconds, with_mpi)

    if clean_workdir:
        builder.clean_workdir = Bool(True)

    launch.launch_process(builder, daemon)
