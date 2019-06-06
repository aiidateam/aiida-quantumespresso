# -*- coding: utf-8 -*-
"""Module with launch utitlies for the CLI."""
from __future__ import absolute_import
import click

from aiida.engine import launch
from .display import echo_calculation_results


def launch_process(process, daemon, **inputs):
    """Launch a process with the given inputs.

    If not sent to the daemon, the results will be displayed after the calculation finishes.

    :param process: the process class
    :param daemon: boolean, if True will submit to the daemon instead of running in current interpreter
    :param inputs: inputs for the process
    """
    if daemon:
        node = launch.submit(process, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(process.__name__, node.pk))
    else:
        click.echo('Running a {}...'.format(process.__name__))
        _, node = launch.run_get_node(process, **inputs)
        echo_calculation_results(node)
