# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PwBandsWorkChain` for testing and demonstration purposes."""
import click

from aiida.cmdline.params import options as options_core
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators

from ...utils import launch
from ...utils import options
from ...utils import validate
from .. import cmd_launch


@cmd_launch.command('pw-bands')
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
@decorators.with_dbenv()
def launch_workflow(
    code, structure, pseudo_family, kpoints_distance, ecutwfc, ecutrho, hubbard_u, hubbard_v, hubbard_file_pk,
    starting_magnetization, smearing, automatic_parallelization, clean_workdir, max_num_machines, max_wallclock_seconds,
    with_mpi, daemon
):
    """Run a `PwBandsWorkChain`."""
    # pylint: disable=too-many-statements
    from aiida.orm import Bool, Float, Dict
    from aiida.plugins import WorkflowFactory
    from qe_tools import CONSTANTS

    from aiida_quantumespresso.utils.resources import get_default_options, get_automatic_parallelization_options

    builder = WorkflowFactory('quantumespresso.pw.bands').get_builder()

    cutoff_wfc, cutoff_rho = pseudo_family.get_recommended_cutoffs(structure=structure)

    parameters = {
        'CONTROL': {
            'calculation': 'relax',
        },
        'SYSTEM': {
            'ecutwfc': ecutwfc or cutoff_wfc / CONSTANTS.ry_to_ev,
            'ecutrho': ecutrho or cutoff_rho / CONSTANTS.ry_to_ev,
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

    pseudos = pseudo_family.get_pseudos(structure=structure)
    parameters = Dict(dict=parameters)

    builder.structure = structure
    builder.relax.base.pw.code = code
    builder.relax.base.pw.parameters = parameters
    builder.relax.base.pw.pseudos = pseudos
    builder.relax.base.kpoints_distance = Float(kpoints_distance)
    builder.relax.meta_convergence = Bool(False)
    builder.scf.pw.code = code
    builder.scf.pw.parameters = parameters
    builder.scf.pw.pseudos = pseudos
    builder.scf.kpoints_distance = Float(kpoints_distance)
    builder.bands.pw.code = code
    builder.bands.pw.parameters = parameters
    builder.bands.pw.pseudos = pseudos

    if hubbard_file:
        builder.relax.base.pw.hubbard_file = hubbard_file
        builder.scf.base.pw.hubbard_file = hubbard_file
        builder.bands.base.pw.hubbard_file = hubbard_file

    if automatic_parallelization:
        auto_para = Dict(dict=get_automatic_parallelization_options(max_num_machines, max_wallclock_seconds))
        builder.relax.base.automatic_parallelization = auto_para
        builder.scf.automatic_parallelization = auto_para
        builder.bands.automatic_parallelization = auto_para
    else:
        metadata_options = get_default_options(max_num_machines, max_wallclock_seconds, with_mpi)
        builder.relax.base.pw.metadata.options = metadata_options
        builder.scf.pw.metadata.options = metadata_options
        builder.bands.pw.metadata.options = metadata_options

    if clean_workdir:
        builder.clean_workdir = Bool(True)

    launch.launch_process(builder, daemon)
