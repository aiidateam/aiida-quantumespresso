# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the `parse_output_base`."""
from pathlib import Path
import pytest


@pytest.fixture
def generate_output(filepath_tests):
    """Return the absolute filepath to the directory containing the file `aiida.out`."""
    fixtures_path = Path(filepath_tests) / 'parsers' / 'fixtures'
    out_file = fixtures_path / 'parse_raw' / 'base' / 'aiida.out'

    with open(out_file) as fout:
        filecontent = fout.read()

    return filecontent


def test_parse_output_base_default(generate_output, data_regression):
    """Test ``parse_output_base`` on the stdout of a ``open_grid.x`` calculation."""
    from aiida_quantumespresso.parsers.parse_raw.base import parse_output_base

    filecontent = generate_output
    parsed_data, logs = parse_output_base(filecontent, codename='OPEN_GRID')

    data_regression.check({'parsed_data': dict(parsed_data), 'logs': dict(logs)})
