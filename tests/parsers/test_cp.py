# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the `CpParser`."""
import pytest

from aiida import orm
from aiida.common import AttributeDict


@pytest.mark.parametrize('version', ['default', '6.6_autopilot', '6.6_verlet', '6.6_cgstep', '6.6_cgsteps'])
def test_cp_default(
    fixture_localhost, generate_calc_job_node, generate_parser, data_regression, generate_structure, version
):
    """Test a default `cp.x` calculation."""
    entry_point_calc_job = 'quantumespresso.cp'
    entry_point_parser = 'quantumespresso.cp'
    if version == 'default':

        def generate_inputs():
            return AttributeDict({
                'structure': generate_structure(structure_id='silicon'),
                'parameters': orm.Dict(dict={}),
            })
    else:

        def generate_inputs():
            return AttributeDict({
                'structure': generate_structure(structure_id='water'),
                'parameters': orm.Dict(dict={}),
            })

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, version, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results

    parameters = results['output_parameters'].get_dict()
    data = {
        'parameters': parameters,
    }

    if version != '6.6_cgstep':
        # In a single cg step we don't have the trajectory output and a message is produced in the log
        assert not orm.Log.objects.get_logs_for(node)
        assert 'output_trajectory' in results
        data['trajectory'] = results['output_trajectory'].attributes

    data_regression.check(data)
