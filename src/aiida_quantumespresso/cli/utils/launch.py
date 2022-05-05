# -*- coding: utf-8 -*-
"""Module with launch utitlies for the CLI."""
import click

from .display import echo_process_results


def launch_process(process, daemon, **inputs):
    """Launch a process with the given inputs.

    If not sent to the daemon, the results will be displayed after the calculation finishes.

    :param process: the process class
    :param daemon: boolean, if True will submit to the daemon instead of running in current interpreter
    :param inputs: inputs for the process
    """
    from aiida.engine import Process, ProcessBuilder, launch

    if isinstance(process, ProcessBuilder):
        process_name = process.process_class.__name__
    elif issubclass(process, Process):
        process_name = process.__name__
    else:
        raise TypeError(f'invalid type for process: {process}')

    if daemon:
        node = launch.submit(process, **inputs)
        click.echo(f'Submitted {process_name}<{node.pk}> to the daemon')
    else:
        if inputs.get('metadata', {}).get('dry_run', False):
            click.echo(f'Running a dry run for {process_name}...')
        else:
            click.echo(f'Running a {process_name}...')
        _, node = launch.run_get_node(process, **inputs)
        echo_process_results(node)
