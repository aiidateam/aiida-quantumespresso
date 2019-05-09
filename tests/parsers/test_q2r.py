# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `Q2rParser`."""
from __future__ import absolute_import

import pytest


@pytest.fixture
def generate_inputs():
    """Return only those inputs that the parser will expect to be there."""
    return {}


def test_q2r_default(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
                     generate_inputs, data_regression):
    """Test a default `q2r.x` calculation."""
    entry_point_calc_job = 'quantumespresso.q2r'
    entry_point_parser = 'quantumespresso.q2r'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default', generate_inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'forceconstants' in results
    data_regression.check(results['forceconstants'].get_content())
