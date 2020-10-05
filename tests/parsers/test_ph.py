# -*- coding: utf-8 -*-
"""Tests for the `PhParser`."""
from aiida import orm


def generate_inputs():
    """Return only those inputs that the parser will expect to be there."""
    return {}


def test_ph_default(fixture_localhost, generate_calc_job_node, generate_parser, data_regression):
    """Test a default `ph.x` calculation."""
    name = 'default'
    entry_point_calc_job = 'quantumespresso.ph'
    entry_point_parser = 'quantumespresso.ph'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())


def test_ph_not_converged(fixture_localhost, generate_calc_job_node, generate_parser, data_regression):
    """Test a `ph.x` calculation where convergence is not reached."""
    name = 'failed_convergence_not_reached'
    entry_point_calc_job = 'quantumespresso.ph'
    entry_point_parser = 'quantumespresso.ph'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_CONVERGENCE_NOT_REACHED.status
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())


def test_ph_out_of_walltime(fixture_localhost, generate_calc_job_node, generate_parser, data_regression):
    """Test a `ph.x` calculation that runs out of walltime."""
    name = 'failed_out_of_walltime'
    entry_point_calc_job = 'quantumespresso.ph'
    entry_point_parser = 'quantumespresso.ph'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUT_OF_WALLTIME.status
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())
