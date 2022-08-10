# -*- coding: utf-8 -*-
# pylint: disable=cyclic-import,reimported,unused-import,wrong-import-position
"""Module with CLI commands for the various work chain implementations."""
from .. import cmd_root


@cmd_root.group('workflow')
def cmd_workflow():
    """Commands to launch and interact with workflows."""


@cmd_workflow.group('launch')
def cmd_launch():
    """Launch workflows."""


# Import the sub commands to register them with the CLI
from .matdyn.base import launch_workflow
from .ph.base import launch_workflow
from .pw.bands import launch_workflow
from .pw.base import launch_workflow
from .pw.relax import launch_workflow
from .q2r.base import launch_workflow
