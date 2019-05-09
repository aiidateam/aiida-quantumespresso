# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `PhParser`."""
from __future__ import absolute_import

import pytest
from aiida import orm
from aiida.common import AttributeDict


@pytest.fixture
def generate_inputs():
    structure = orm.StructureData()
    parameters = {'INPUTPH': {}}
    qpoints = orm.KpointsData()
    qpoints.set_cell_from_structure(structure)
    qpoints.set_kpoints_mesh([2, 2, 2])

    return AttributeDict({
        'parent_folder': orm.FolderData(),
        'qpoints': qpoints,
        'parameters': orm.Dict(dict=parameters),
    })


def test_ph_default(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
                    generate_inputs, data_regression):
    """Test a default `ph.x` calculation.

    The output is created by running a dead simple SCF calculation for a silicon structure.
    This test should test the standard parsing of the stdout content and XML file stored in the standard results node.
    """
    entry_point_calc_job = 'quantumespresso.ph'
    entry_point_parser = 'quantumespresso.ph'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default', generate_inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())
