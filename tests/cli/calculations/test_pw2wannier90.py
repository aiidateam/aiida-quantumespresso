# -*- coding: utf-8 -*-
"""Tests for the ``calculation launch pw2wannier90`` command."""
import io

from aiida.orm import SinglefileData

from aiida_quantumespresso.cli.calculations.pw2wannier90 import launch_calculation


def test_command_base(run_cli_process_launch_command, fixture_code, generate_calc_job_node):
    """Test invoking the calculation launch command with only required inputs."""
    code = fixture_code('quantumespresso.pw2wannier90').store()
    calculation = generate_calc_job_node('quantumespresso.pw', test_name='default').store()
    nnkp_file = SinglefileData(io.BytesIO(b'content')).store()
    options = ['-X', code.full_label, '-P', calculation.outputs.remote_folder.pk, '-S', nnkp_file.pk]
    run_cli_process_launch_command(launch_calculation, options=options)
