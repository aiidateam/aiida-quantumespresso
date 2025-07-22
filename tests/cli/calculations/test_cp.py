# -*- coding: utf-8 -*-
"""Tests for the ``calculation launch cp`` command."""
from aiida_quantumespresso.cli.calculations.cp import launch_calculation


def test_command_base(run_cli_process_launch_command, fixture_code, pseudo_family):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.cp').store()
    options = ['-X', code.full_label, '-F', pseudo_family.label]
    run_cli_process_launch_command(launch_calculation, options=options)
