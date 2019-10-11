# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `PwParser`."""
from __future__ import absolute_import

import pytest

from aiida import orm
from aiida.common import AttributeDict
import numpy as np


@pytest.fixture
def generate_neb_structures():
    """Return 2 `StructureData` objects that can be used as first and last images for a NEB calculation."""

    def _generate_neb_structures(element='H'):
        """Return 2 `StructureData` objects that can be used as first and last images for a NEB calculation."""
        from aiida.orm import StructureData
        from qe_tools.constants import bohr_to_ang
        import numpy as np

        cell = np.array([[12, 0, 0], [0, 5, 0], [0, 0, 5]]) * bohr_to_ang

        structure_1 = StructureData(cell=cell)
        structure_1.append_atom(position=np.array([-4.56670009, 0., 0.]) * bohr_to_ang, symbols=element)
        structure_1.append_atom(position=(0., 0., 0.), symbols=element)
        structure_1.append_atom(position=np.array([1.55776676, 0., 0.]) * bohr_to_ang, symbols=element)

        structure_2 = StructureData(cell=cell)
        structure_2.append_atom(position=np.array([-1.55776676, 0., 0.]) * bohr_to_ang, symbols=element)
        structure_2.append_atom(position=(0., 0., 0.), symbols=element)
        structure_2.append_atom(position=np.array([4.56670009, 0., 0.]) * bohr_to_ang, symbols=element)

        return structure_1, structure_2

    return _generate_neb_structures


@pytest.fixture
def generate_inputs(generate_neb_structures):

    def _generate_inputs(ci_scheme='auto', parser_options=None):
        """Return only those inputs that the parser will expect to be there."""
        first_structure, last_structure = generate_neb_structures()

        if ci_scheme in ['auto', 'no-CI']:
            num_images = 3
        elif ci_scheme == 'manual':
            num_images = 4
        else:
            raise ValueError('Unexpected ci_scheme')

        settings = {
            'fixed_coords':
                [[False, True, True],
                 [True, True, True],
                 [False, True, True]],
            'parser_options': parser_options,
        }  # yapf: disable
        if num_images == 'manual':
            settings['CLIMBING_IMAGES'] = [num_images // 2 + 1]

        neb_parameters = {
            'PATH': {
                'nstep_path': 20,
                'ds': 2.,
                'opt_scheme': 'broyden',
                'num_of_images': num_images,
                'k_max': 0.3,
                'k_min': 0.2,
                'CI_scheme': ci_scheme,
                'path_thr': 0.05,
            },
        }

        pw_parameters = {
            'CONTROL': {
                'calculation': 'relax'
            },
            'SYSTEM': {
                'ecutwfc': 20,
                'ecutrho': 100,
                'occupations': 'smearing',
                'degauss': 0.003,
                'nspin': 2,
                'starting_magnetization': 0.5,
            },
            'ELECTRONS': {
                'conv_thr': 1e-8,
                'mixing_beta': 0.3,
            }
        }

        kpoints = orm.KpointsData()
        kpoints.set_kpoints_mesh([2, 2, 2])

        inputs = {
            'parameters': orm.Dict(dict=neb_parameters),
            'pw': {
                'parameters': orm.Dict(dict=pw_parameters),
                'kpoints': kpoints,
            },
            'first_structure': first_structure,
            'last_structure': last_structure,
            'settings': orm.Dict(dict=settings)
        }
        return AttributeDict(inputs)

    return _generate_inputs


def test_neb_default(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
                     generate_inputs, data_regression, num_regression):
    """Test a NEB calculation with symmetric images and automatic climbing image."""
    name = 'default'
    entry_point_calc_job = 'quantumespresso.neb'
    entry_point_parser = 'quantumespresso.neb'

    inputs = generate_inputs(ci_scheme='auto', parser_options=None)
    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_mep' in results
    assert 'output_trajectory' in results
    assert 'iteration_array' not in results

    data_regression.check({
        'parameters': results['output_parameters'].get_dict()
    })

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
    """Test a NEB calculation with symmetric images and automatic climbing image."""
    name = 'default'
    entry_point_calc_job = 'quantumespresso.neb'
    entry_point_parser = 'quantumespresso.neb'

    inputs = generate_inputs(ci_scheme='auto', parser_options={'all_iterations': True})
    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_mep' in results
    assert 'output_trajectory' in results
    assert 'iteration_array' in results

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
    """Test a NEB calculation with symmetric images and automatic climbing image."""
    name = 'default'
    entry_point_calc_job = 'quantumespresso.neb'
    entry_point_parser = 'quantumespresso.neb'

    inputs = generate_inputs(ci_scheme='auto', parser_options={'include_deprecated_v2_keys': True})
    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_mep' in results
    assert 'output_trajectory' in results
    assert 'iteration_array' not in results

    for img in [1, 2, 3]:
        pw_params = results['output_parameters']['pw_output_image_{}'.format(img)]
        assert pw_params['fixed_occupations'] == False
        assert pw_params['smearing_method'] == True
        assert pw_params['tetrahedron_method'] == False
