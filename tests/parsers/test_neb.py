# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `NebParser`."""
from __future__ import absolute_import

import pytest

from aiida import orm
from aiida.common import AttributeDict
import numpy as np


@pytest.fixture
def generate_inputs():

    def _generate_inputs(parser_options=None):
        """Return only those inputs that the parser will expect to be there."""
        inputs = {
            'parameters': orm.Dict(dict={'PATH': {'num_of_images': 3}}),
            'pw': {
                'parameters': orm.Dict(dict={}),
            },
            'settings': orm.Dict(dict={'parser_options': parser_options})
        }
        return AttributeDict(inputs)

    return _generate_inputs


def test_neb_default(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
        generate_inputs, data_regression, num_regression):
    """Test a NEB calculation with symmetric images and automatic climbing image."""
    name = 'default'
    entry_point_calc_job = 'quantumespresso.neb'
    entry_point_parser = 'quantumespresso.neb'

    inputs = generate_inputs(parser_options=None)
    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
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

    num_data_dict = {}

    num_data_dict.update({
        arr: results['output_mep'].get_array(arr).flatten() for arr in ['mep', 'interpolated_mep']
    })
    num_data_dict.update({
        arr: results['output_trajectory'].get_array(arr).flatten() for arr in ['cells', 'positions']
    })

    # Convert all arrays to floats, to get around this change that disallows diffent-sized arrays for non-float types:
    # https://github.com/ESSS/pytest-regressions/pull/18
    for key, val in num_data_dict.items():
        if not (np.issubdtype(val.dtype, np.floating) or np.issubdtype(val.dtype, np.complexfloating)):
            num_data_dict[key] = val.astype(np.float64)

    num_regression.check(num_data_dict, default_tolerance=dict(atol=0, rtol=1e-18))


def test_neb_all_iterations(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
        generate_inputs, data_regression, num_regression):
    """Test a NEB calculation with the parser option `all_iterations=True`."""
    name = 'default'
    entry_point_calc_job = 'quantumespresso.neb'
    entry_point_parser = 'quantumespresso.neb'

    inputs = generate_inputs(parser_options={'all_iterations': True})
    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_mep' in results
    assert 'output_trajectory' in results
    assert 'iteration_array' in results

    data = {
        'iteration_array': results['iteration_array'].attributes,
    }
    data_regression.check(data)

    num_data_dict = {
        arr: results['iteration_array'].get_array(arr).flatten() for arr in results['iteration_array'].get_arraynames()
    }

    # Convert all arrays to floats, to get around this change that disallows diffent-sized arrays for non-float types:
    # https://github.com/ESSS/pytest-regressions/pull/18
    for key, val in num_data_dict.items():
        if not (np.issubdtype(val.dtype, np.floating) or np.issubdtype(val.dtype, np.complexfloating)):
            num_data_dict[key] = val.astype(np.float64)

    num_regression.check(num_data_dict, default_tolerance=dict(atol=0, rtol=1e-18))


def test_neb_deprecated_keys(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
        generate_inputs, data_regression, num_regression):
    """Test a NEB calculation with the parser option `include_deprecated_v2_keys=True`."""
    name = 'default'
    entry_point_calc_job = 'quantumespresso.neb'
    entry_point_parser = 'quantumespresso.neb'

    inputs = generate_inputs(parser_options={'include_deprecated_v2_keys': True})
    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_mep' in results
    assert 'output_trajectory' in results
    assert 'iteration_array' not in results

    for key, dictionary in results['output_parameters'].get_dict().items():
        if key.startswith('pw_output_image'):
            assert dictionary['fixed_occupations'] is False
            assert dictionary['smearing_method'] is True
            assert dictionary['tetrahedron_method'] is False
