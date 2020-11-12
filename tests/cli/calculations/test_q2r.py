# -*- coding: utf-8 -*-
"""Tests for the ``calculation launch q2r`` command."""
from aiida_quantumespresso.cli.calculations.q2r import launch_calculation


def test_command_base(run_cli_process_launch_command, fixture_code, generate_calc_job_node):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.q2r').store()
    calculation = generate_calc_job_node('quantumespresso.ph', test_name='default').store()
    options = ['-X', code.full_label, '-C', calculation.pk]
    run_cli_process_launch_command(launch_calculation, options=options)
