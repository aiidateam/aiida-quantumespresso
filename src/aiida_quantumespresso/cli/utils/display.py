# -*- coding: utf-8 -*-
"""Module with display utitlies for the CLI."""
import os

import click


def echo_process_results(node):
    """Display a formatted table of the outputs registered for the given process node.

    :param node: the `ProcessNode` of a terminated process
    """
    from aiida.common.links import LinkType

    class_name = node.process_class.__name__
    outputs = node.base.links.get_outgoing(link_type=(LinkType.CREATE, LinkType.RETURN)).all()

    if hasattr(node, 'dry_run_info'):
        # It is a dry-run: get the information and print it
        rel_path = os.path.relpath(node.dry_run_info['folder'])
        click.echo(f"-> Files created in folder '{rel_path}'")
        click.echo(f"-> Submission script filename: '{node.dry_run_info['script_filename']}'")
        return

    if node.is_finished and node.exit_message:
        state = f'{node.process_state.value} [{node.exit_status}] `{node.exit_message}`'
    elif node.is_finished:
        state = f'{node.process_state.value} [{node.exit_status}]'
    else:
        state = node.process_state.value

    click.echo(f'{class_name}<{node.pk}> terminated with state: {state}')

    if not outputs:
        click.echo(f'{class_name}<{node.pk}> registered no outputs')
        return

    click.echo(f"\n{'Output link':25s} Node pk and type")
    click.echo(f"{'-' * 60}")

    for triple in sorted(outputs, key=lambda triple: triple.link_label):
        click.echo(f'{triple.link_label:25s} {triple.node.__class__.__name__}<{triple.node.pk}> ')
