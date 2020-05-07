# -*- coding: utf-8 -*-
"""Module with display utitlies for the CLI."""
from __future__ import absolute_import

import os

import click


def echo_process_results(node):
    """Display a formatted table of the outputs registered for the given process node.

    :param node: the `ProcessNode` of a terminated process
    """
    from aiida.common.links import LinkType

    class_name = node.process_class.__name__
    outputs = node.get_outgoing(link_type=(LinkType.CREATE, LinkType.RETURN)).all()

    if hasattr(node, 'dry_run_info'):
        # It is a dry-run: get the information and print it
        rel_path = os.path.relpath(node.dry_run_info['folder'])
        click.echo("-> Files created in folder '{}'".format(rel_path))
        click.echo("-> Submission script filename: '{}'".format(node.dry_run_info['script_filename']))
        return

    if node.is_finished and node.exit_message:
        state = '{} [{}] `{}`'.format(node.process_state.value, node.exit_status, node.exit_message)
    elif node.is_finished:
        state = '{} [{}]'.format(node.process_state.value, node.exit_status)
    else:
        state = node.process_state.value

    click.echo('{}<{}> terminated with state: {}'.format(class_name, node.pk, state))

    if not outputs:
        click.echo('{}<{}> registered no outputs'.format(class_name, node.pk))
        return

    click.echo('\n{link:25s} {node}'.format(link='Output link', node='Node pk and type'))
    click.echo('{s}'.format(s='-' * 60))

    for triple in sorted(outputs, key=lambda triple: triple.link_label):
        click.echo('{:25s} {}<{}> '.format(triple.link_label, triple.node.__class__.__name__, triple.node.pk))
