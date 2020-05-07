# -*- coding: utf-8 -*-
"""Test the automatic 'invalidates_cache' attribute for exit codes."""
import inspect
import pytest

from aiida.engine import CalcJob
from aiida.plugins import CalculationFactory
from aiida.plugins.entry_point import get_entry_point_names

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
    entry_point = CalculationFactory(entry_point_name)

    if not inspect.isclass(entry_point) or not issubclass(entry_point, CalcJob):
        return

    overrides = EXPLICIT_OVERRIDES.get(entry_point_name, [])

    for exit_code in entry_point.exit_codes.values():
        if exit_code.status not in overrides:
            if exit_code.status < 400:
                assert exit_code.invalidates_cache
            else:
                assert not exit_code.invalidates_cache
