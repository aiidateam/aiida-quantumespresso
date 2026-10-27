"""Module for the command line interface."""

import click
from aiida.cmdline.groups import VerdiCommandGroup
from aiida.cmdline.params import options, types

from .setup.codes import setup_codes_cmd


@click.group(
    'aiida-quantumespresso',
    cls=VerdiCommandGroup,
    context_settings={'help_option_names': ['-h', '--help']},
)
@options.PROFILE(type=types.ProfileParamType(load_profile=True), expose_value=False)
def cmd_root():
    """CLI for the `aiida-quantumespresso` plugin."""


@cmd_root.group('setup')
def cmd_setup():
    """Easy setup of codes, computers, and other essentials."""


cmd_setup.command('codes')(setup_codes_cmd)

from .calculations import cmd_calculation
from .data import cmd_structure
from .workflows import cmd_workflow
