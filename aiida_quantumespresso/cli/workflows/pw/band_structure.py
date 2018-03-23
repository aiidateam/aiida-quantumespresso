# -*- coding: utf-8 -*-
import click
from aiida.utils.cli import command
from aiida.utils.cli import options


@command()
@options.code(callback_kwargs={'entry_point': 'quantumespresso.pw'})
@options.structure()
@options.pseudo_family()
@options.daemon()
@click.option(
    '-z', '--protocol', type=click.Choice(['standard']), default='standard', show_default=True,
    help='the protocol to use for the workflow'
)
def launch(
    code, structure, pseudo_family, daemon, protocol):
    """
    Run the PwBandStructureWorkChain for a given input structure 
    to compute the band structure for the relaxed structure
    """
    from aiida.orm.data.base import Str
    from aiida.orm.utils import WorkflowFactory
    from aiida.work.launch import run, submit

    PwBandStructureWorkChain = WorkflowFactory('quantumespresso.pw.band_structure')

    inputs = {
        'code': code,
        'structure': structure,
        'pseudo_family': Str(pseudo_family),
        'protocol': Str(protocol),
    }

    if daemon:
        workchain = submit(PwBandStructureWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(PwBandStructureWorkChain.__name__, workchain.pk))
    else:
        run(PwBandStructureWorkChain, **inputs)
