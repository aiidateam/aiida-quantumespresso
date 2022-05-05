# -*- coding: utf-8 -*-
# pylint: disable=cyclic-import,reimported,unused-import,wrong-import-position
"""Module with CLI commands for the various calculation job implementations."""
from .. import cmd_root


@cmd_root.group('calculation')
def cmd_calculation():
    """Commands to launch and interact with calculations."""


@cmd_calculation.group('launch')
def cmd_launch():
    """Launch calculations."""


# Import the sub commands to register them with the CLI
from .cp import launch_calculation
from .dos import launch_calculation
from .epw import launch_calculation
from .matdyn import launch_calculation
from .neb import launch_calculation
from .ph import launch_calculation
from .pp import launch_calculation
from .projwfc import launch_calculation
from .pw2wannier90 import launch_calculation
from .pw import launch_calculation
from .q2r import launch_calculation
