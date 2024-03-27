# -*- coding: utf-8 -*-
"""Tests for the ``data structure import`` command."""
from pathlib import Path
import re

import pytest

from aiida_quantumespresso.cli.data.structure import cmd_import


@pytest.mark.usefixtures('aiida_profile')
def test_command_import(run_cli_command, filepath_tests):
    """Test invoking the calculation launch command with only required inputs."""
    filepath = str(Path(filepath_tests, 'cli', 'fixtures', 'pw.in').absolute())
    options = [filepath]
    result = run_cli_command(cmd_import, options=options)
    match = r'.*parsed and stored StructureData<(\d+)> with formula (.+?)'
    assert re.match(match, ' '.join(result.output_lines))

    options = [filepath, '--dry-run']
    result = run_cli_command(cmd_import, options=options)
    match = r'.*parsed structure with formula (.+?)'
    assert re.match(match, ' '.join(result.output_lines))
