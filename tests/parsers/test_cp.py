# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the `CpParser`."""
import pytest

from aiida import orm
from aiida.common import AttributeDict


@pytest.fixture
def generate_inputs(generate_structure):
    """Return only those inputs that the parser will expect to be there."""
    return AttributeDict({
        'structure': generate_structure(),
        'parameters': orm.Dict(dict={}),
    })


def test_cp_default(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression):
    """Test a default `cp.x` calculation."""
    entry_point_calc_job = 'quantumespresso.cp'
    entry_point_parser = 'quantumespresso.cp'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default', generate_inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'parameters': results['output_parameters'].get_dict(),
        'trajectory': results['output_trajectory'].attributes
    })
