# -*- coding: utf-8 -*-
"""Tests for the ``calculation launch epw`` command."""
from aiida_quantumespresso.cli.calculations.epw import launch_calculation


def test_command_base(run_cli_process_launch_command, fixture_code, generate_calc_job_node):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.epw').store()
    parent_pw = generate_calc_job_node('quantumespresso.pw', test_name='default').store()
    parent_ph = generate_calc_job_node('quantumespresso.ph', test_name='default').store()
    options = ['-X', code.full_label, '--pw-nscf-parent', parent_pw.pk, '--ph-parent', parent_ph.pk]
    run_cli_process_launch_command(launch_calculation, options=options)
