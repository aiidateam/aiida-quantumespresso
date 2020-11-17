# -*- coding: utf-8 -*-
"""Tests for the ``workflow launch pw-base`` command."""
from aiida_quantumespresso.cli.workflows.pw.base import launch_workflow


def test_command_base(run_cli_process_launch_command, fixture_code, sssp):
    """Test invoking the workflow launch command with only required inputs."""
    code = fixture_code('quantumespresso.pw').store()
    options = ['-X', code.full_label, '-F', sssp.label]
    run_cli_process_launch_command(launch_workflow, options=options)
