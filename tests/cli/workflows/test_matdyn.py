# -*- coding: utf-8 -*-
"""Tests for the ``calculation launch matdyn`` command."""
from aiida_quantumespresso.cli.calculations.matdyn import launch_calculation


def test_command_base(run_cli_process_launch_command, fixture_code, generate_force_constants_data):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.matdyn').store()
    force_constants = generate_force_constants_data.store()
    options = ['-X', code.full_label, '-D', force_constants.pk]
    run_cli_process_launch_command(launch_calculation, options=options)
