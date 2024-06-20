# -*- coding: utf-8 -*-
"""Tests for the `PhParser`."""
from aiida import orm
import pytest


def generate_inputs():
    """Return only those inputs that the parser will expect to be there."""
    return {'parameters': orm.Dict({'INPUTPH': {}})}


@pytest.mark.parametrize('test_name', ['default', 'single_qpoint', 'no_modes_printed'])
def test_ph_default(test_name, fixture_localhost, generate_calc_job_node, generate_parser, data_regression):
    """Default tests for the `ph.x` parser."""
    name = test_name
    entry_point_calc_job = 'quantumespresso.ph'
    entry_point_parser = 'quantumespresso.ph'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not [log for log in orm.Log.collection.get_logs_for(node) if log.levelname == 'ERROR']
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


def test_ph_initialization(fixture_localhost, generate_calc_job_node, generate_parser, data_regression):
    """Test a `ph.x` calculation performed with `start_irr` and `last_irr` set to 0."""
    name = 'initialization'
    entry_point_calc_job = 'quantumespresso.ph'
    entry_point_parser = 'quantumespresso.ph'

    inputs = {'parameters': orm.Dict({'INPUTPH': {'start_irr': 0, 'last_irr': 0}})}

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())


def test_ph_initialization_failed(fixture_localhost, generate_calc_job_node, generate_parser):
    """Test a failed `ph.x` calculation performed with `start_irr` and `last_irr` set to 0."""
    name = 'initialization_failed'
    entry_point_calc_job = 'quantumespresso.ph'
    entry_point_parser = 'quantumespresso.ph'

    inputs = {'parameters': orm.Dict({'INPUTPH': {'start_irr': 0, 'last_irr': 0}})}

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE.status


def test_ph_failed_computing_cholesky(fixture_localhost, generate_calc_job_node, generate_parser):
    """Test the parsing of a calculation that failed during cholesky factorization.

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = 'failed_computing_cholesky'
    entry_point_calc_job = 'quantumespresso.ph'
    entry_point_parser = 'quantumespresso.ph'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_COMPUTING_CHOLESKY.status


def test_ph_failed_incompatible_fft(fixture_localhost, generate_calc_job_node, generate_parser):
    """Test the parsing of a calculation that failed finding an incompatible FFT grid."""
    name = 'failed_incompatible_fft'
    entry_point_calc_job = 'quantumespresso.ph'
    entry_point_parser = 'quantumespresso.ph'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_INCOMPATIBLE_FFT_GRID.status


def test_ph_failed_wrong_representation(fixture_localhost, generate_calc_job_node, generate_parser):
    """Test the parsing of a calculation that failed finding a wrong representation."""
    name = 'failed_wrong_representation'
    entry_point_calc_job = 'quantumespresso.ph'
    entry_point_parser = 'quantumespresso.ph'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_WRONG_REPRESENTATION.status
