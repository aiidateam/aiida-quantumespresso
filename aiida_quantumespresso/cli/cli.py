# -*- coding: utf-8 -*-
"""Module for the command line interface."""
from __future__ import absolute_import
import click


@click.group('aiida-quantumespresso', context_settings={'help_option_names': ['-h', '--help']})
def root():
    """CLI for the `aiida-quantumespresso` plugin."""


@root.group('calculation')
def calculation():
    """Commands to launch and interact with calculations."""


@calculation.group('launch')
def calculation_launch():
    """Launch calculations."""


@root.group('workflow')
def workflow():
    """Commands to launch and interact with workflows."""


@workflow.group('launch')
def workflow_launch():
    """Launch workflows."""
