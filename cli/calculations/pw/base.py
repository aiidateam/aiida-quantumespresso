# -*- coding: utf-8 -*-
import click
from cli.utils import options
from cli.utils.click import command


@command()
@options.code()
@options.structure()
@options.pseudo_family()
@options.kpoint_mesh()
@options.max_num_machines()
@options.max_wallclock_seconds()
@click.option(
    '-z', '--calculation-mode', 'mode', type=click.Choice(['scf', 'vc-relax']), default='scf', show_default=True,
    help='select the calculation mode'
)
def launch(code, structure, pseudo_family, kpoints, max_num_machines, max_wallclock_seconds, mode):
    """
    Run a PwCalculation for a given input structure
    """
    from aiida.orm import load_node
    from aiida.orm.data.parameter import ParameterData
    from aiida.orm.data.upf import get_pseudos_from_structure
    from aiida.orm.utils import CalculationFactory
    from aiida.work.run import run
    from aiida_quantumespresso.utils.resources import get_default_options

    PwCalculation = CalculationFactory('quantumespresso.pw')

    parameters = {
        'CONTROL': {
            'calculation': mode,
        },
        'SYSTEM': {
            'ecutwfc': 30.,
            'ecutrho': 240.,
        },
    }

    inputs = {
        'code': code,
        'structure': structure,
        'pseudo': get_pseudos_from_structure(structure, pseudo_family),
        'kpoints': kpoints,
        'parameters': ParameterData(dict=parameters),
        'settings': ParameterData(dict={}),
        '_options': get_default_options(max_num_machines, max_wallclock_seconds),
    }

    click.echo('Running a pw.x calculation in the {} mode... '.format(mode))

    process = PwCalculation.process()
    results, pk = run(process, _return_pid=True, **inputs)
    calculation = load_node(pk)

    click.echo('PwCalculation<{}> terminated with state: {}'.format(pk, calculation.get_state()))
    click.echo('\n{link:25s} {node}'.format(link='Output link', node='Node pk and type'))
    click.echo('{s}'.format(s='-'*60))
    for link, node in sorted(calculation.get_outputs(also_labels=True)):
        click.echo('{:25s} <{}> {}'.format(link, node.pk, node.__class__.__name__))