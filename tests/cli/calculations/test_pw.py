# -*- coding: utf-8 -*-
"""Tests for the ``calculation launch pw`` command."""
import re

import pytest

from aiida_quantumespresso.cli.calculations.pw import launch_calculation


def test_command_base(run_cli_process_launch_command, fixture_code, pseudo_family):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.pw').store()
    options = ['-X', code.full_label, '-F', pseudo_family.label]
    run_cli_process_launch_command(launch_calculation, options=options)


# yapf: disable
@pytest.mark.parametrize(('cmd_options', 'match'), (
    (
        ['--hubbard-u', 'Mg'],
        ".*Option '--hubbard-u' requires 2 arguments.*"),
    (
        ['--hubbard-u', 'Mg', '5.0'],
        '.*kinds in the specified Hubbard U is not a strict subset of the structure kinds.*'),
    (
        ['--hubbard-file', '1000000'],
        '.*no SinglefileData found with ID*'
    ),
))
# yapf: enable
def test_invalid_hubbard_parameters(run_cli_process_launch_command, fixture_code, pseudo_family, cmd_options, match):
    """Test invoking the calculation launch command with invalid Hubbard inputs."""
    code = fixture_code('quantumespresso.pw').store()
    options = ['-X', code.full_label, '-F', pseudo_family.label] + cmd_options
    result = run_cli_process_launch_command(launch_calculation, options=options, raises=ValueError)
    assert re.match(match, ' '.join(result.output_lines))


def test_valid_hubbard_parameters(run_cli_process_launch_command, fixture_code, pseudo_family):
    """Test invoking the calculation launch command with valid Hubbard inputs."""
    import io

    from aiida.orm import SinglefileData

    code = fixture_code('quantumespresso.pw').store()

    options = ['-X', code.full_label, '-F', pseudo_family.label, '--hubbard-u', 'Si', '5.0']
    run_cli_process_launch_command(launch_calculation, options=options)

    content_original = 'for sure some correct Hubbard parameters'
    filepk = SinglefileData(io.StringIO(content_original)).store().pk

    options = ['-X', code.full_label, '-F', pseudo_family.label, '--hubbard-file', filepk]
    run_cli_process_launch_command(launch_calculation, options=options)
