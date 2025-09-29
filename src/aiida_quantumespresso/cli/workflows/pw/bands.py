"""Command line scripts to launch a `PwBandsWorkChain` for testing and demonstration purposes."""

import click
from aiida.cmdline.params import options as options_core
from aiida.cmdline.params import types
from aiida.cmdline.utils import decorators

from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain

from ...utils import launch, options, validate
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
@options.CLEAN_WORKDIR()
@options.MAX_NUM_MACHINES()
@options.MAX_WALLCLOCK_SECONDS()
@options.WITH_MPI()
@options.DAEMON()
@decorators.with_dbenv()
def launch_workflow(
    code,
    structure,
    pseudo_family,
    kpoints_distance,
    ecutwfc,
    ecutrho,
    hubbard_u,
    hubbard_v,
    hubbard_file,
    starting_magnetization,
    smearing,
    clean_workdir,
    max_num_machines,
    max_wallclock_seconds,
    with_mpi,
    daemon,
):
    """Run a `PwBandsWorkChain`."""
    from aiida_quantumespresso.utils.resources import get_default_options

    cutoff_wfc, cutoff_rho = pseudo_family.get_recommended_cutoffs(structure=structure, unit='Ry')

    parameters = {
        'SYSTEM': {
            'ecutwfc': ecutwfc or cutoff_wfc,
            'ecutrho': ecutrho or cutoff_rho,
        },
    }

    try:
        validate.validate_hubbard_parameters(structure, parameters, hubbard_u, hubbard_v, hubbard_file)
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

    base_overrides = {
        'clean_workdir': clean_workdir,
        'pseudo_family': pseudo_family.label,
        'kpoints_distance': kpoints_distance,
        'pw': {
            'parameters': parameters,
            'metadata': {'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi)},
            'hubbard_file': hubbard_file,
        },
    }
    overrides = {'relax': base_overrides, 'scf': base_overrides, 'bands': base_overrides}
    builder = PwBandsWorkChain.get_builder_from_protocol(code, structure, overrides=overrides)

    launch.launch_process(builder, daemon)
