# -*- coding: utf-8 -*-
"""Test the automatic 'invalidates_cache' attribute for exit codes."""

from __future__ import absolute_import

from distutils.version import StrictVersion  # pylint: disable=import-error,no-name-in-module

import pytest

import aiida
from aiida.plugins import CalculationFactory
from aiida.plugins.entry_point import get_entry_point_names

if StrictVersion(aiida.__version__) < StrictVersion('1.1.0'):
    pytest.skip("The 'invalidates_cache' feature is only available on AiiDA 1.1+", allow_module_level=True)

QE_CALC_ENTRY_POINT_NAMES = [
    ep_name for ep_name in get_entry_point_names(group='aiida.calculations') if ep_name.startswith('quantumespresso')
]

# When explicitly overriden 'invalidates_cache' are added, add entries
# EXPLICIT_OVERRIDES = {<entry_point_name>: [<status_integer>, ...]}

EXPLICIT_OVERRIDES = {}


@pytest.mark.parametrize('entry_point_name', QE_CALC_ENTRY_POINT_NAMES)
def test_exit_code_invalidates_cache(entry_point_name):
    """Test automatic 'invalidates_cache' attribute of exit codes.

    Test that the 'invalidates_cache' attribute of exit codes is automatically
    set according to the status integer.
    """
    calc_class = CalculationFactory(entry_point_name)
    overrides = EXPLICIT_OVERRIDES.get(entry_point_name, [])
    for exit_code in calc_class.exit_codes.values():
        if exit_code.status not in overrides:
            if exit_code.status < 400:
                assert exit_code.invalidates_cache
            else:
                assert not exit_code.invalidates_cache
