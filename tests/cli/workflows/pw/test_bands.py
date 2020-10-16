# -*- coding: utf-8 -*-
"""Tests for the ``workflow launch pw-bands`` command."""
import pytest

from aiida_quantumespresso.cli.workflows.pw.bands import launch_workflow


@pytest.mark.usefixtures('clear_database_before_test')
def test_command_base(run_cli_process_launch_command, fixture_code, generate_upf_family, generate_structure):
    """Test invoking the workflow launch command with only required inputs."""
    code = fixture_code('quantumespresso.pw').store()
    family = generate_upf_family()
    options = ['-X', code.full_label, '-p', family.label, '-s', generate_structure().store().pk]
    run_cli_process_launch_command(launch_workflow, options=options)
