# -*- coding: utf-8 -*-
"""Tests for the `NebParser`."""
import numpy as np

from aiida import orm
from aiida.common import AttributeDict


def generate_inputs(parser_options=None):
    """Return only those inputs that the parser will expect to be there."""
    inputs = {
        'parameters': orm.Dict(dict={'PATH': {
            'num_of_images': 3
        }}),
        'pw': {
            'parameters': orm.Dict()
        },
        'settings': orm.Dict(dict={'parser_options': parser_options})
    }
    return AttributeDict(inputs)


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


def test_neb_default(fixture_localhost, generate_calc_job_node, generate_parser, data_regression, num_regression):
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
    assert not orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_mep' in results
    assert 'output_trajectory' in results
    assert 'iteration_array' not in results

    data = {
        'parameters': results['output_parameters'].get_dict(),
        'output_mep': results['output_mep'].attributes,
        'output_trajectory': results['output_trajectory'].attributes,
    }
    data_regression.check(data)

    data = build_num_regression_dictionary([results['output_mep']], [['mep', 'interpolated_mep']])
    num_regression.check(data, default_tolerance=dict(atol=0, rtol=1e-18))


def test_neb_all_iterations(
    fixture_localhost, generate_calc_job_node, generate_parser, data_regression, num_regression
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
    assert not orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_mep' in results
    assert 'output_trajectory' in results
    assert 'iteration_array' in results

    data = {'iteration_array': results['iteration_array'].attributes}
    data_regression.check(data)

    data = build_num_regression_dictionary([results['iteration_array']], [results['iteration_array'].get_arraynames()])
    num_regression.check(data, default_tolerance=dict(atol=0, rtol=1e-18))
