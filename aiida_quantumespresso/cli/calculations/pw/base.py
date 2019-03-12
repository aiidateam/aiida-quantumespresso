# -*- coding: utf-8 -*-
import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from aiida_quantumespresso.cli.utils import options as options_qe
from aiida_quantumespresso.cli.utils import validate


@click.command()
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.pw'))
@options_qe.STRUCTURE(required=True)
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
@click.option(
    '-z', '--calculation-mode', 'mode', type=click.Choice(['scf', 'vc-relax']), default='scf', show_default=True,
    help='select the calculation mode'
)
@decorators.with_dbenv()
def launch(
    code, structure, pseudo_family, kpoints_mesh, ecutwfc, ecutrho, hubbard_u, hubbard_v, hubbard_file_pk,
    starting_magnetization, smearing, max_num_machines, max_wallclock_seconds, with_mpi, daemon, mode):
    """Run a PwCalculation."""
    from aiida.orm import Dict
    from aiida.orm.nodes.data.upf import get_pseudos_from_structure
    from aiida.engine import launch
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    PwCalculation = CalculationFactory('quantumespresso.pw')

    parameters = {
        'CONTROL': {
            'calculation': mode,
        },
        'SYSTEM': {
            'ecutwfc': ecutwfc,
            'ecutrho': ecutrho,
        },
    }

    try:
        hubbard_file = validate.validate_hubbard_parameters(structure, parameters, hubbard_u, hubbard_v, hubbard_file_pk)
    except ValueError as exception:
        raise click.BadParameter(exception.message)

    try:
        validate.validate_starting_magnetization(structure, parameters, starting_magnetization)
    except ValueError as exception:
        raise click.BadParameter(exception.message)

    try:
        validate.validate_smearing(parameters, smearing)
    except ValueError as exception:
        raise click.BadParameter(exception.message)

    inputs = {
        'code': code,
        'structure': structure,
        'pseudos': get_pseudos_from_structure(structure, pseudo_family),
        'kpoints': kpoints_mesh,
        'parameters': Dict(dict=parameters),
        'settings': Dict(dict={}),
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    if hubbard_file:
        inputs['hubbard_file'] = hubbard_file

    if daemon:
        calculation = launch.submit(PwCalculation, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(PwCalculation.__name__, calculation.pk))
    else:
        click.echo('Running a pw.x calculation in the {} mode... '.format(mode))
        results, calculation = launch.run_get_node(PwCalculation, **inputs)
        click.echo('PwCalculation<{}> terminated with state: {}'.format(calculation.pk, calculation.get_state()))
        click.echo('\n{link:25s} {node}'.format(link='Output link', node='Node pk and type'))
        click.echo('{s}'.format(s='-' * 60))
        for triple in sorted(calculation.get_outgoing().all()):
            click.echo('{:25s} {}<{}> '.format(triple.link_label, triple.node.__class__.__name__, triple.node.pk))
