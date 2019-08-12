# -*- coding: utf-8 -*-
# pylint: disable=wrong-import-position
"""Module with CLI commands for various data types."""
from ..cli import root


@root.group('data')
def cmd_data():
    """Commands to create and inspect data nodes."""


from .structure import cmd_structure
