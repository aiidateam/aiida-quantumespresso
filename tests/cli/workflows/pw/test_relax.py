# -*- coding: utf-8 -*-
"""Tests for the ``workflow launch pw-relax`` command."""
from aiida_quantumespresso.cli.workflows.pw.relax import launch_workflow


def test_command_base(run_cli_process_launch_command, fixture_code, generate_structure, sssp):
    """Test invoking the workflow launch command with only required inputs."""
    code = fixture_code('quantumespresso.pw').store()
    options = ['-X', code.full_label, '-F', sssp.label, '-S', generate_structure().store().pk]
    run_cli_process_launch_command(launch_workflow, options=options)
