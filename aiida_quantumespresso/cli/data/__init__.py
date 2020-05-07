# -*- coding: utf-8 -*-
# pylint: disable=cyclic-import,unused-import,wrong-import-position
"""Module with CLI commands for various data types."""
from .. import cmd_root


@cmd_root.group('data')
def cmd_data():
    """Commands to create and inspect data nodes."""


# Import the sub commands to register them with the CLI
from .structure import cmd_structure
