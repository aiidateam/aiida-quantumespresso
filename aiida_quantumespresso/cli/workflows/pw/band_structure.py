# -*- coding: utf-8 -*-
import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from aiida_quantumespresso.cli.utils import options as options_qe


@click.command()
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.pw'))
@options_qe.STRUCTURE(required=True)
@options_qe.DAEMON()
@click.option(
    '-z', '--protocol', type=click.Choice(['standard']), default='standard', show_default=True,
    help='the protocol to use for the workflow'
)
@decorators.with_dbenv()
def launch(
    code, structure, pseudo_family, daemon, protocol):
    """Run a PwBandStructureWorkChain."""
    from aiida.plugins import DataFactory
    from aiida.plugins import WorkflowFactory
    from aiida.work import launch

    PwBandStructureWorkChain = WorkflowFactory('quantumespresso.pw.band_structure')
    Dict = DataFactory('dict')

    inputs = {
        'code': code,
        'structure': structure,
        'protocol': Dict(dict={
            'name': 'theos-ht-1.0',
        }),
    }

    if daemon:
        workchain = launch.submit(PwBandStructureWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(PwBandStructureWorkChain.__name__, workchain.pk))
    else:
        launch.run(PwBandStructureWorkChain, **inputs)
