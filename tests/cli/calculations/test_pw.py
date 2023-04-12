# -*- coding: utf-8 -*-
"""Tests for the ``calculation launch pw`` command."""
import re

import pytest

from aiida_quantumespresso.cli.calculations.pw import launch_calculation


def test_command_base(run_cli_process_launch_command, fixture_code, sssp):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.pw').store()
    options = ['-X', code.full_label, '-F', sssp.label]
    run_cli_process_launch_command(launch_calculation, options=options)


@pytest.mark.parametrize(('cmd_options', 'match'), (
    (['--hubbard-u', 'Mg'], ".*Option '--hubbard-u' requires 2 arguments.*"),
    (['--hubbard-u', 'Mg', '5.0'
      ], '.*kinds in the specified Hubbard U is not a strict subset of the structure kinds.*'),
    (['--hubbard-file', '1000000'], '.*no SinglefileData found with ID*'),
))
def test_invalid_hubbard_parameters(run_cli_process_launch_command, fixture_code, sssp, cmd_options, match):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.pw').store()
    options = ['-X', code.full_label, '-F', sssp.label] + cmd_options
    result = run_cli_process_launch_command(launch_calculation, options=options, raises=ValueError)
    assert re.match(match, ' '.join(result.output_lines))


@pytest.mark.usefixtures('aiida_profile')
def test_valid_hubbard_parameters(run_cli_process_launch_command, fixture_code, sssp):
    """Test invoking the calculation launch command with only required inputs."""
    import pathlib
    import tempfile

    from aiida.orm import SinglefileData

    code = fixture_code('quantumespresso.pw').store()

    options = ['-X', code.full_label, '-F', sssp.label, '--hubbard-u', 'Si', '5.0']
    run_cli_process_launch_command(launch_calculation, options=options)

    content_original = 'a very creative content'

    with tempfile.NamedTemporaryFile(mode='w+') as handle:
        filepath = pathlib.Path(handle.name).resolve()
        handle.write(content_original)
        handle.flush()
        filepk = SinglefileData(file=filepath).store().pk

    options = ['-X', code.full_label, '-F', sssp.label, '--hubbard-file', filepk]
    run_cli_process_launch_command(launch_calculation, options=options)
