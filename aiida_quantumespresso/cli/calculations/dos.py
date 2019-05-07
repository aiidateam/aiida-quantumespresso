# -*- coding: utf-8 -*-
"""Command line scripts to launch a `PwCalculation` for testing and demonstration purposes."""
from __future__ import absolute_import
import click

from aiida.orm import RemoteData
from aiida.cmdline.params import options, types
from aiida.cmdline.utils import decorators

from aiida_quantumespresso.cli.utils import options as options_qe


@click.command()
@options.CODE(required=True, type=types.CodeParamType(entry_point='quantumespresso.dos'))
@options.CALCULATION(required=True)
@options_qe.MAX_NUM_MACHINES()
@options_qe.MAX_WALLCLOCK_SECONDS()
@options_qe.WITH_MPI()
@options_qe.DAEMON()
@decorators.with_dbenv()
def cli(code, calculation, max_num_machines, max_wallclock_seconds, with_mpi, daemon):
    """Run a DosCalculation."""
    from aiida.engine import launch
    from aiida.plugins import CalculationFactory
    from aiida_quantumespresso.utils.resources import get_default_options

    DosCalculation = CalculationFactory('quantumespresso.dos')  # pylint: disable=invalid-name

    # Check that the parent calculation node comes from quantumespresso.pw.
    # I cannot move this check into the option declaration, because CalcJobNode is not subclassed by the specific
    # calculation plugins (only Process is), and there is no feature yet to filter by the associated process_type.
    expected_process_type = 'aiida.calculations:quantumespresso.pw'
    if calculation.process_type != expected_process_type:
        raise click.BadParameter('The input calculation node has a process_type: {}; should be {}'.format(
            calculation.process_type, expected_process_type))
    parent_folder = calculation.get_outgoing(node_class=RemoteData, link_label_filter='remote_folder').one().node

    inputs = {
        'code': code,
        'parent_folder': parent_folder,
        'metadata': {
            'options': get_default_options(max_num_machines, max_wallclock_seconds, with_mpi),
        }
    }

    if daemon:
        new_calc = launch.submit(DosCalculation, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(DosCalculation.__name__, new_calc.pk))
    else:
        click.echo('Running a dos.x calculation from parent {}<{}>... '.format(calculation.__class__.__name__,
                                                                               calculation.pk))
        _, new_calc = launch.run_get_node(DosCalculation, **inputs)
        click.echo('DosCalculation<{}> terminated with state: {}'.format(new_calc.pk, new_calc.process_state))
        click.echo('\n{link:25s} {node}'.format(link='Output link', node='Node pk and type'))
        click.echo('{s}'.format(s='-' * 60))
        for triple in sorted(new_calc.get_outgoing().all(), key=lambda triple: triple.link_label):
            click.echo('{:25s} {}<{}> '.format(triple.link_label, triple.node.__class__.__name__, triple.node.pk))
