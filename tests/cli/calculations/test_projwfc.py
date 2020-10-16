# -*- coding: utf-8 -*-
"""Tests for the ``calculation launch projwfc`` command."""
import pytest

from aiida_quantumespresso.cli.calculations.projwfc import launch_calculation


@pytest.mark.usefixtures('clear_database_before_test')
def test_command_base(run_cli_process_launch_command, fixture_code, generate_calc_job_node):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.projwfc').store()
    calculation = generate_calc_job_node('quantumespresso.pw', test_name='default').store()
    options = ['-X', code.full_label, '-C', calculation.pk]
    run_cli_process_launch_command(launch_calculation, options=options)
