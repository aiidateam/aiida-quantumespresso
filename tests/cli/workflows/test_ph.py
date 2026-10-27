"""Tests for the ``workflows launch ph-base`` command."""

from aiida_quantumespresso.cli.workflows.ph.base import launch_workflow


def test_command_base(run_cli_process_launch_command, fixture_code, generate_calc_job_node):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.ph').store()
    calculation = generate_calc_job_node('quantumespresso.pw', test_name='default').store()
    options = ['-X', code.full_label, '-C', calculation.pk]
    run_cli_process_launch_command(launch_workflow, options=options)
