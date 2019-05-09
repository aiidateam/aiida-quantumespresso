# -*- coding: utf-8 -*-
"""Command line scripts to launch a `MatdynCalculation` for testing and demonstration purposes."""
from __future__ import absolute_import
import click

from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from aiida_quantumespresso.cli.utils import options as options_qe


@click.command()
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.matdyn'))
@options.DATUM(
    required=True,
    type=types.DataParamType(sub_classes=('aiida.data:quantumespresso.forceconstants',)),
    help='A ForceconstantsData node produced by a `Q2rCalculation`')
@options_qe.KPOINTS_MESH(default=[1, 1, 1])
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def cli(code, datum, kpoints_mesh, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run a MatdynCalculation."""
    from aiida.engine import launch
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    MatdynCalculation = CalculationFactory('quantumespresso.matdyn')  # pylint: disable=invalid-name

    inputs = {
        'code': code,
        'kpoints': kpoints_mesh,
        'parent_folder': datum,
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    if daemon:
        node = launch.submit(MatdynCalculation, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(MatdynCalculation.__name__, node.pk))
    else:
        click.echo('Running a matdyn.x calculation from {}<{}>... '.format(datum.__class__.__name__, datum.pk))
        _, node = launch.run_get_node(MatdynCalculation, **inputs)
        click.echo('MatdynCalculation<{}> terminated with state: {}'.format(node.pk, node.process_state))
        click.echo('\n{link:25s} {node}'.format(link='Output link', node='Node pk and type'))
        click.echo('{s}'.format(s='-' * 60))
        for triple in sorted(node.get_outgoing().all(), key=lambda triple: triple.link_label):
            click.echo('{:25s} {}<{}> '.format(triple.link_label, triple.node.__class__.__name__, triple.node.pk))
