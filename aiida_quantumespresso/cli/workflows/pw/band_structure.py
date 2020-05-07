# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PwBandStructureWorkChain` for testing and demonstration purposes."""
import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from ...utils import launch
from ...utils import options as options_qe
from .. import cmd_launch


@cmd_launch.command('pw-band-structure')
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.pw'))
@options_qe.STRUCTURE(required=True)
@options_qe.DAEMON()
@click.option(
    '-z',
    '--protocol',
    type=click.Choice(['theos-ht-1.0']),
    default='theos-ht-1.0',
    show_default=True,
    help='The protocol to use for the workflow.'
)
@decorators.with_dbenv()
def launch_workflow(code, structure, daemon, protocol):
    """Run a `PwBandStructureWorkChain`."""
    from aiida import orm
    from aiida.plugins import WorkflowFactory

    inputs = {
        'code': code,
        'structure': structure,
        'protocol': orm.Dict(dict={'name': protocol}),
    }

    launch.launch_process(WorkflowFactory('quantumespresso.pw.band_structure'), daemon, **inputs)
