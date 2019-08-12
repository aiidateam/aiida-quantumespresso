# -*- coding: utf-8 -*-
"""Module with launch utitlies for the CLI."""
from __future__ import absolute_import

import click

from .display import echo_process_results


def launch_process(process, daemon, **inputs):
    """Launch a process with the given inputs.

    If not sent to the daemon, the results will be displayed after the calculation finishes.

    :param process: the process class
    :param daemon: boolean, if True will submit to the daemon instead of running in current interpreter
    :param inputs: inputs for the process
    """
    from aiida.engine import launch, Process, ProcessBuilder

    if isinstance(process, ProcessBuilder):
        process_name = process.process_class.__name__
    elif issubclass(process, Process):
        process_name = process.__name__
    else:
        raise TypeError('invalid type for process: {}'.format(process))

    if daemon:
        node = launch.submit(process, **inputs)
        click.echo('Submitted {}<{}> to the daemon'.format(process_name, node.pk))
    else:
        if inputs.get('metadata', {}).get('dry_run', False):
            click.echo('Running a dry run for {}...'.format(process_name))
        else:
            click.echo('Running a {}...'.format(process_name))
        _, node = launch.run_get_node(process, **inputs)
        echo_process_results(node)
