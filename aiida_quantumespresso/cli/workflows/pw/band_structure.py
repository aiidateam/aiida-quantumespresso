# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PwBandStructureWorkChain` for testing and demonstration purposes."""
import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from aiida_quantumespresso.cli.utils import options as options_qe


@click.command()
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.pw'))
@options_qe.STRUCTURE(required=True)
@options_qe.DAEMON()
@click.option(
    '-z',
    '--protocol',
    type=click.Choice(['theos-ht-1.0']),
    default='theos-ht-1.0',
    show_default=True,
    help='the protocol to use for the workflow')
@decorators.with_dbenv()
def cli(code, structure, daemon, protocol):
    """Run a `PwBandStructureWorkChain`."""
    from aiida.engine import launch
    from aiida.plugins import DataFactory, WorkflowFactory

    PwBandStructureWorkChain = WorkflowFactory('quantumespresso.pw.band_structure')  # pylint: disable=invalid-name
    Dict = DataFactory('dict')  # pylint: disable=invalid-name

    inputs = {
        'code': code,
        'structure': structure,
        'protocol': Dict(dict={'name': protocol}),
    }

    if daemon:
        workchain = launch.submit(PwBandStructureWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(PwBandStructureWorkChain.__name__, workchain.pk))
    else:
        launch.run(PwBandStructureWorkChain, **inputs)
