# -*- coding: utf-8 -*-
import click
from aiida.utils.cli import command
from aiida.utils.cli import options
from aiida_quantumespresso.utils.cli import options as options_qe


@command()
@options.code()
@options.structure()
@options.pseudo_family()
@options.kpoint_mesh()
@options.max_num_machines()
@options.max_wallclock_seconds()
@options.daemon()
@options_qe.ecutwfc()
@options_qe.ecutrho()
@options_qe.hubbard_u()
@click.option(
    '-z', '--calculation-mode', 'mode', type=click.Choice(['scf', 'vc-relax']), default='scf', show_default=True,
    help='Select the calculation mode'
)
def launch(code, structure, pseudo_family, kpoints, max_num_machines, max_wallclock_seconds, daemon, mode, hubbard_u,
    ecutwfc, ecutrho):
    """
    Run a PwCalculation for a given input structure
    """
    from aiida.orm import load_node
    from aiida.orm.data.parameter import ParameterData
    from aiida.orm.data.upf import get_pseudos_from_structure
    from aiida.orm.utils import CalculationFactory
    from aiida.work.launch import run_get_node, submit
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

    if hubbard_u:
        structure_kinds = structure.get_kind_names()
        hubbard_u_kinds = [kind for kind, value in hubbard_u]

        if not set(hubbard_u_kinds).issubset(structure_kinds):
            raise click.BadParameter('the kinds in the specified starting Hubbard U values is not a strict subset of the kinds in the structure')
            return

        parameters['SYSTEM']['hubbard_u'] = {}
        parameters['SYSTEM']['lda_plus_u'] = True

        for kind, value in hubbard_u:
            parameters['SYSTEM']['hubbard_u'][kind] = value

    inputs = {
        'code': code,
        'structure': structure,
        'pseudo': get_pseudos_from_structure(structure, pseudo_family),
        'kpoints': kpoints,
        'parameters': ParameterData(dict=parameters),
        'settings': ParameterData(dict={}),
        'options': get_default_options(max_num_machines, max_wallclock_seconds),
    }


    process = PwCalculation.process()

    if daemon:
        calculation = submit(process, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(PwCalculation.__name__, calculation.pk))
    else:
        click.echo('Running a pw.x calculation in the {} mode... '.format(mode))
        results, calculation = run_get_node(process, **inputs)
        click.echo('PwCalculation<{}> terminated with state: {}'.format(calculation.pk, calculation.get_state()))
        click.echo('\n{link:25s} {node}'.format(link='Output link', node='Node pk and type'))
        click.echo('{s}'.format(s='-'*60))
        for link, node in sorted(calculation.get_outputs(also_labels=True)):
            click.echo('{:25s} <{}> {}'.format(link, node.pk, node.__class__.__name__))
