# -*- coding: utf-8 -*-
"""Tests for the ``calculation launch neb`` command."""
from aiida_quantumespresso.cli.calculations.neb import launch_calculation


def test_command_base(run_cli_process_launch_command, fixture_code, sssp, generate_structure):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.neb').store()
    structures = [generate_structure().store().pk, generate_structure().store().pk]
    options = ['-X', code.full_label, '-F', sssp.label, '-s'] + structures
    run_cli_process_launch_command(launch_calculation, options=options)
