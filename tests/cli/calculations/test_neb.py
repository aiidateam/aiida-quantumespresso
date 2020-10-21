# -*- coding: utf-8 -*-
"""Tests for the ``calculation launch neb`` command."""
import pytest

from aiida_quantumespresso.cli.calculations.neb import launch_calculation


@pytest.mark.usefixtures('clear_database_before_test')
def test_command_base(run_cli_process_launch_command, fixture_code, generate_upf_family, generate_structure):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.neb').store()
    family = generate_upf_family()
    structures = [generate_structure().store().pk, generate_structure().store().pk]
    options = ['-X', code.full_label, '-p', family.label, '-s'] + structures
    run_cli_process_launch_command(launch_calculation, options=options)
