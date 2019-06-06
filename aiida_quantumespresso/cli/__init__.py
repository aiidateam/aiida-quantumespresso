# -*- coding: utf-8 -*-
# pylint: disable=wrong-import-position,wildcard-import
"""Module for the command line interface."""
from __future__ import absolute_import
import click_completion

# Activate the completion of parameter types provided by the click_completion package
click_completion.init()

# Import to populate the `verdi` sub commands
from .cli import calculation, calculation_launch, workflow, workflow_launch
from .calculations import *
from .workflows import *
