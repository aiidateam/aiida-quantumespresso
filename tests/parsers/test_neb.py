# -*- coding: utf-8 -*-
# pylint: disable=invalid-name,redefined-outer-name, too-many-lines
"""Tests for the `NebParser`."""
from aiida import orm
from aiida.common import AttributeDict
import numpy as np
import pytest

from aiida_quantumespresso.calculations.neb import NebCalculation


@pytest.fixture
def generate_inputs(generate_trajectory):
    """Generate the inputs for the NEB parser."""

    def _generate_inputs(parser_options=None):
        """Return only those inputs that the parser will expect to be there."""
        inputs = {
            'images': generate_trajectory(),
            'parameters': orm.Dict({'PATH': {
                'num_of_images': 3
            }}),
            'pw': {
                'parameters': orm.Dict()
            },
            'settings': orm.Dict({'parser_options': parser_options})
        }
        return AttributeDict(inputs)

    return _generate_inputs


def build_num_regression_dictionary(arrays, array_names):
    """Build a dictionary that can be passed to `num_regression`.

    :param arrays: a list of `ArrayData` nodes
    :param array_names: a list of array name lists, should have same length as `arrays`
    :return: dictionary
    """
    if len(arrays) != len(array_names):
        raise ValueError('length of `arrays` and `array_names` should be equal.')

    result = {}

    for index, array in enumerate(arrays):
        for name in array_names[index]:
            result[name] = array.get_array(name).flatten()

    # Convert all arrays to floats, to get around this change that disallows diffent-sized arrays for non-float types:
    # https://github.com/ESSS/pytest-regressions/pull/18
    for key, val in result.items():
        if not (np.issubdtype(val.dtype, np.floating) or np.issubdtype(val.dtype, np.complexfloating)):  # pylint: disable=no-member
            result[key] = val.astype(np.float64)

    return result


def test_neb_default(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression, num_regression
):
    """Test a NEB calculation with symmetric images and automatic climbing image."""
    name = 'default'
    entry_point_calc_job = 'quantumespresso.neb'
    entry_point_parser = 'quantumespresso.neb'

    inputs = generate_inputs(parser_options=None)
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not [log for log in orm.Log.collection.get_logs_for(node) if log.levelname == 'ERROR']
    assert not [
        log for log in orm.Log.collection.get_logs_for(node)
        if 'DEPRECATED: symmetry with ibrav=0, use correct ibrav instead' not in log.message
    ]
    assert 'output_parameters' in results
    assert 'output_mep' in results
    assert 'output_trajectory' in results
    assert 'iteration_array' not in results

    data = {
        'parameters': results['output_parameters'].get_dict(),
        'output_mep': results['output_mep'].base.attributes.all,
        'output_trajectory': results['output_trajectory'].base.attributes.all,
    }
    data_regression.check(data)

    data = build_num_regression_dictionary([results['output_mep']], [['mep', 'interpolated_mep']])
    num_regression.check(data, default_tolerance=dict(atol=0, rtol=1e-18))


def test_neb_all_iterations(
    fixture_localhost, generate_calc_job_node, generate_parser, data_regression, num_regression, generate_inputs
):
    """Test a NEB calculation with the parser option `all_iterations=True`."""
    name = 'default'
    entry_point_calc_job = 'quantumespresso.neb'
    entry_point_parser = 'quantumespresso.neb'

    inputs = generate_inputs(parser_options={'all_iterations': True})
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not [log for log in orm.Log.collection.get_logs_for(node) if log.levelname == 'ERROR']
    assert 'output_parameters' in results
    assert 'output_mep' in results
    assert 'output_trajectory' in results
    assert 'iteration_array' in results

    data = {'iteration_array': results['iteration_array'].base.attributes.all}
    data_regression.check(data)

    data = build_num_regression_dictionary([results['iteration_array']], [results['iteration_array'].get_arraynames()])
    num_regression.check(data, default_tolerance=dict(atol=0, rtol=1e-18))


@pytest.mark.parametrize(
    'filename, exception', [
        ('failed_computing_cholesky', 'ERROR_COMPUTING_CHOLESKY'),
        ('failed_too_many_bands_not_converged', 'ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED'),
        ('failed_s_matrix_not_positive_definite', 'ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE'),
        ('failed_zhegvd_failed', 'ERROR_ZHEGVD_FAILED'),
        ('failed_qr_failed', 'ERROR_QR_FAILED'),
        ('failed_eigenvectors_convergence', 'ERROR_EIGENVECTOR_CONVERGENCE'),
        ('failed_broyden_factorization', 'ERROR_BROYDEN_FACTORIZATION'),
        ('failed_dexx_negative', 'ERROR_DEXX_IS_NEGATIVE'),
    ]
)
def test_failed_diagonalization(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename, exception
):
    """Test the parsing of a calculation that failed with diagonalization exceptions.

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """

    entry_point_calc_job = 'quantumespresso.neb'
    entry_point_parser = 'quantumespresso.neb'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, filename, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.get(exception).status


@pytest.mark.parametrize(
    'test_case, expected_exit_code', (
        ('default', None),
        ('failed_interrupted', NebCalculation.exit_codes.ERROR_NEB_INTERRUPTED_PARTIAL_TRAJECTORY),
    )
)
def test_failed_interrupted_scheduler(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, test_case, expected_exit_code
):
    """Test that an exit code set by the scheduler is not overridden unless a more specific error is parsed.

    The test is run twice, once for the ``default`` test case and once for ``failed_interrupted``, which correspond to a
    successful run and a run that got interrupted (usually due to scheduler killing the job). Before calling the parser
    the ``ERROR_SCHEDULER_OUT_OF_WALLTIME`` is set on the node. In the case of the ``default`` test case, this should be
    ignored and the parser should return ``ExitCode(0)``. For the interrupted case, the exit code of the scheduler
    should not be overridden so the parser should return ``ERROR_NEB_INTERRUPTED_PARTIAL_TRAJECTORY``.
    """
    entry_point_calc_job = 'quantumespresso.neb'
    entry_point_parser = 'quantumespresso.neb'

    # Generate the node and set an exit status as if it would have been set by the scheduler parser
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, test_case, generate_inputs())
    node.set_exit_status(NebCalculation.exit_codes.ERROR_SCHEDULER_OUT_OF_WALLTIME.status)

    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    if expected_exit_code is None:
        assert calcfunction.is_finished_ok, (calcfunction.exit_status, calcfunction.exception)
    else:
        assert not calcfunction.is_finished_ok, calcfunction.exception
        assert calcfunction.exit_status == expected_exit_code.status


def test_failed_cycle_exceeded_nstep(
    fixture_localhost,
    generate_calc_job_node,
    generate_parser,
    generate_inputs,
):
    """Test the parsing of a calculation that stopped due to exceeding the maximum number of steps."""
    name = 'failed_cycle_exceeded_nstep'
    entry_point_calc_job = 'quantumespresso.neb'
    entry_point_parser = 'quantumespresso.neb'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_NEB_CYCLE_EXCEEDED_NSTEP.status
