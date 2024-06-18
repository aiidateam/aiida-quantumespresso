# -*- coding: utf-8 -*-
"""Tests for the `BandsParser`."""
from aiida import orm


def generate_inputs():
    """Return only those inputs that the parser will expect to be there."""
    return {}


def test_bands_default(fixture_localhost, generate_calc_job_node, generate_parser, data_regression):
    """Test a default `bands.x` calculation."""
    entry_point_calc_job = 'quantumespresso.bands'
    entry_point_parser = 'quantumespresso.bands'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default', generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.collection.get_logs_for(node)
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())
