# -*- coding: utf-8 -*-
import click
from aiida_quantumespresso.utils.click import command
from aiida_quantumespresso.utils.click import options


@command()
@options.code()
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
    from aiida.work.run import run, submit

    PwBandStructureWorkChain = WorkflowFactory('quantumespresso.pw.band_structure')

    inputs = {
        'code': code,
        'structure': structure,
        'pseudo_family': Str(pseudo_family),
        'protocol': Str(protocol),
    }

    if daemon:
        workchain = submit(PwBandStructureWorkChain, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(PwBandStructureWorkChain.__name__, workchain.pid))
    else:
        run(PwBandStructureWorkChain, **inputs)
