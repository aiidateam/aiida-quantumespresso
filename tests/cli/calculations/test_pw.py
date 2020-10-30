# -*- coding: utf-8 -*-
"""Tests for the ``calculation launch pw`` command."""
import pytest

from aiida_quantumespresso.cli.calculations.pw import launch_calculation


@pytest.mark.usefixtures('clear_database_before_test')
def test_command_base(run_cli_process_launch_command, fixture_code, generate_upf_family):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.pw').store()
    family = generate_upf_family()
    options = ['-X', code.full_label, '-F', family.label]
    run_cli_process_launch_command(launch_calculation, options=options)
