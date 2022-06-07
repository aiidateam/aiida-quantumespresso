# -*- coding: utf-8 -*-
"""Tests for the :class:`aiida_quantumespresso.calculations.importers.pw.PwCalculationImporter` class."""
from pathlib import Path

from aiida.engine import run
from aiida.orm import RemoteData

from aiida_quantumespresso.calculations.importers.pw import PwCalculationImporter
from aiida_quantumespresso.calculations.pw import PwCalculation


def test_default(filepath_tests, fixture_code, aiida_localhost):
    """Test importing a typical completed ``pw.x`` calculation."""
    aiida_localhost.configure()
    filepath_remote = Path(filepath_tests) / 'calculations' / 'importers' / 'test_pw' / 'test_default'
    remote_data = RemoteData(str(filepath_remote), computer=aiida_localhost)
    input_file_name = 'aiida.in'
    code = fixture_code('quantumespresso.pw')
    pseudo_folder_path = filepath_remote / 'pseudo'
    inputs = PwCalculationImporter().parse_remote_data(
        remote_data, input_file_name, code, pseudo_folder_path=str(pseudo_folder_path)
    )
    results, node = run.get_node(PwCalculation, **inputs)
    assert node.is_finished_ok
    assert node.is_imported
    assert set(results.keys()) == {'output_band', 'output_trajectory', 'output_parameters', 'retrieved'}
