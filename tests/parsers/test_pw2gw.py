# -*- coding: utf-8 -*-
# pylint: disable=invalid-name,redefined-outer-name
"""Tests for the `Pw2gwParser`."""
from aiida import orm


def test_pw2gw_default(fixture_localhost, generate_parser, generate_calc_job_node, data_regression, num_regression):
    """Test a normal pw2gw.x output."""
    name = 'default'
    entry_point_calc_job = 'quantumespresso.pw2gw'
    entry_point_parser = 'quantumespresso.pw2gw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name)

    parser = generate_parser(entry_point_parser)

    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)

    data_regression.check({'output_parameters': results['output_parameters'].get_dict()},
                          basename='test_pw2gw_default_data')

    num_regression.check(dict(results['eps'].get_iterarrays()), basename='test_pw2gw_default_eps')


def test_pw2gw_failed_missing_output(fixture_localhost, generate_parser, generate_calc_job_node):
    """Test a pw2gw.x output where file are missing."""
    name = 'failed_missing_output'
    entry_point_calc_job = 'quantumespresso.pw2gw'
    entry_point_parser = 'quantumespresso.pw2gw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name)

    parser = generate_parser(entry_point_parser)

    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status, node.process_class.exit_codes.ERROR_OUTPUT_FILES.status
    assert orm.Log.objects.get_logs_for(node)


def test_pw2gw_failed_missing_stdout(fixture_localhost, generate_parser, generate_calc_job_node):
    """Test a pw2gw.x output where file are missing."""
    name = 'failed_missing_stdout'
    entry_point_calc_job = 'quantumespresso.pw2gw'
    entry_point_parser = 'quantumespresso.pw2gw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name)

    parser = generate_parser(entry_point_parser)

    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status, node.process_class.exit_codes.ERROR_OUTPUT_STDOUT_MISSING.status
    assert orm.Log.objects.get_logs_for(node)


def test_pw2gw_failed_corrupted_file(fixture_localhost, generate_parser, generate_calc_job_node):
    """Test a pw2gw.x output where file are corrupted."""
    name = 'failed_corrupted_file'
    entry_point_calc_job = 'quantumespresso.pw2gw'
    entry_point_parser = 'quantumespresso.pw2gw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name)

    parser = generate_parser(entry_point_parser)

    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status, node.process_class.exit_codes.ERROR_OUTPUT_FILES.status
    assert orm.Log.objects.get_logs_for(node)
