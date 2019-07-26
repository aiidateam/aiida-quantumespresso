# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `PwParser`."""
from __future__ import absolute_import

import pytest
from aiida import orm
from aiida.common import AttributeDict


@pytest.fixture
def generate_neb_structures():
    """Return 2 `StructureData` objects that can be used as first and last images for a NEB calculation."""

    def _generate_neb_structures(element='H'):
        """Return 2 `StructureData` objects that can be used as first and last images for a NEB calculation."""
        from aiida.orm import StructureData
        from aiida_quantumespresso.parsers.constants import bohr_to_ang
        import numpy as np

        cell = np.array([[12, 0, 0], [0, 5, 0], [0, 0, 5]]) * bohr_to_ang

        structure_1 = StructureData(cell=cell)
        structure_1.append_atom(position=np.array([-4.56670009, 0., 0.])*bohr_to_ang, symbols=element)
        structure_1.append_atom(position=(0., 0., 0.), symbols=element)
        structure_1.append_atom(position=np.array([1.55776676, 0., 0.])*bohr_to_ang, symbols=element)

        structure_2 = StructureData(cell=cell)
        structure_2.append_atom(position=np.array([-1.55776676, 0., 0.])*bohr_to_ang, symbols=element)
        structure_2.append_atom(position=(0., 0., 0.), symbols=element)
        structure_2.append_atom(position=np.array([4.56670009, 0., 0.])*bohr_to_ang, symbols=element)

        return structure_1, structure_2

    return _generate_neb_structures


@pytest.fixture
def generate_inputs(generate_neb_structures):
    """Return only those inputs that the parser will expect to be there."""
    first_structure, last_structure = generate_neb_structures()

    parameters = {
        'PATH': {
            'ds': 2.0,
            'k_max': 0.3,
            'k_min': 0.2,
            'path_thr': 0.1, #0.05
            'CI_scheme': 'auto', #'manual'
            'nstep_path': 20,
            'opt_scheme': 'broyden',
            'num_of_images': 7, #8
        }
    }

    pw_parameters = {
        'SYSTEM': {
            'nspin': 2,
            'degauss': 0.003,
            'ecutrho': 100.0,
            'ecutwfc': 20.0,
            'occupations': 'smearing',
            'starting_magnetization': 0.5
        },
        'CONTROL': {
            'calculation': 'relax'
        },
        'ELECTRONS': {
            'conv_thr': 1e-08, 'mixing_beta': 0.3
        }
    }

    settings = {
        'fixed_coords':
            [[False, True, True],
             [True, True, True],
             [False, True, True]],
        #'CLIMBING_IMAGES': [5],
    }

    kpoints = orm.KpointsData()
    #kpoints.set_cell_from_structure(structure)
    kpoints.set_kpoints_mesh([2,2,2])

    return AttributeDict({
        'parameters': orm.Dict(dict=parameters),
        'pw': {
            'parameters': orm.Dict(dict=pw_parameters),
            'kpoints': kpoints,
        },
        'first_structure': first_structure,
        'last_structure': last_structure,
        'settings': orm.Dict(dict=settings)
    })


def test_neb_h2h(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
                 generate_inputs, data_regression):
    """Test a default `pw.x` calculation.

    The output is created by running a dead simple SCF calculation for a silicon structure.
    This test should test the standard parsing of the stdout content and XML file stored in the standard results node.
    """
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    # inputs = generate_inputs
    # first_structure, last_structure = generate_neb_structures()
    # inputs.first_structure = first_structure
    # inputs.last_structure = last_structure
    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default', generate_inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())
