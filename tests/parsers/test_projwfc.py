# -*- coding: utf-8 -*-
# pylint: disable=unused-argument,redefined-outer-name
"""Tests for the `ProjwfcParser`."""
from __future__ import absolute_import

import os
import pytest

from aiida import orm
from aiida.common import AttributeDict, LinkType

from aiida_quantumespresso.parsers.parse_raw.projwfc import parse_lowdin_charges


@pytest.mark.parametrize('test_file', ('default', 'spin'))
def test_parse_lowdin_charges(test_file, parser_fixture_path, data_regression):
    path = os.path.join(parser_fixture_path, 'projwfc', test_file, 'aiida.out')
    with open(path) as handle:
        data, spill_param = parse_lowdin_charges(handle.read().splitlines())
    data['spill'] = spill_param
    data_regression.check(data)


@pytest.fixture
def projwfc_inputs(generate_calc_job_node, fixture_computer_localhost, generate_structure, generate_kpoints_mesh):
    """Create the required inputs for the ``ProjwfcCalculation``."""
    parent_calcjob = generate_calc_job_node(
        'quantumespresso.pw',
        fixture_computer_localhost,
        'default',
        inputs={
            'structure': generate_structure('Si'),
            'kpoints': generate_kpoints_mesh(4)
        }
    )
    params = orm.Dict(dict={'number_of_spin_components': 1})
    params.add_incoming(parent_calcjob, link_type=LinkType.CREATE, link_label='output_parameters')
    params.store()
    inputs = {
        'parent_folder': parent_calcjob.outputs.remote_folder,
    }

    return AttributeDict(inputs)


def test_projwfc_default(
    fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser, projwfc_inputs,
    data_regression
):
    """Test ``ProjwfcParser`` on the results of a simple ``projwfc.x`` calculation."""
    entry_point_calc_job = 'quantumespresso.projwfc'
    entry_point_parser = 'quantumespresso.projwfc'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default', projwfc_inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message

    for link_name in ['output_parameters', 'Dos', 'bands', 'projections']:
        assert link_name in results, list(results.keys())

    data_regression.check({
        'Dos': results['Dos'].attributes,
        'bands': results['bands'].attributes,
        'projections':
        {k: v for k, v in results['projections'].attributes.items() if k not in ['reference_bandsdata_uuid']}
    })
