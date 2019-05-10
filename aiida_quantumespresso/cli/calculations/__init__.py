# -*- coding: utf-8 -*-
# pylint: disable=reimported
"""Module with CLI commands for the various calculation job implementations."""
from .cp import launch_calculation
from .dos import launch_calculation
from .matdyn import launch_calculation
from .ph import launch_calculation
from .pp import launch_calculation
from .pw import launch_calculation
from .pw2wannier90 import launch_calculation
from .projwfc import launch_calculation
from .q2r import launch_calculation
