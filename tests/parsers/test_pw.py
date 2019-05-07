# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `PwParser`."""
from __future__ import absolute_import


def test_pw_default(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser, 
                    generate_inputs, data_regression):
    """Test a default `pw.x` calculation.

    The output is created by running a dead simple SCF calculation for a silicon structure.
    This test should test the standard parsing of the stdout content and XML file stored in the standard results node.
    """
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    input_nodes = generate_inputs(entry_point_calc_job)

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default', input_nodes)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())
