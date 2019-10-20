# -*- coding: utf-8 -*-
# pylint: disable=unused-argument,redefined-outer-name
"""Tests for the `ProjwfcParser`."""
from __future__ import absolute_import

import os
import pytest

from aiida import orm
from aiida.common import AttributeDict, LinkType
from aiida.plugins import CalculationFactory

from aiida_quantumespresso.parsers.parse_raw.projwfc import parse_lowdin_charges


@pytest.mark.parametrize('test_file', ('default', 'spin'))
def test_parse_lowdin_charges(test_file, parser_fixture_path, data_regression):
    """Test parsing of lowdin charges from projwfc.x, for non-spin/spin cases."""
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
    parameters = {'PROJWFC': {'emin': -1, 'emax': 1, 'DeltaE': 0.01, 'ngauss': 0, 'degauss': 0.01}}
    inputs = {
        'parameters': orm.Dict(dict=parameters),
        'parent_folder': parent_calcjob.outputs.remote_folder,
    }

    return AttributeDict(inputs)


def test_pw_link_spec():
    """Test the ``PwCalculation`` input/output link names are as required for ``ProjwfcParser``.

    ``ProjwfcParser`` relies on extracting data from the parent ``PwCalculation``.
    This test safeguards against changes in the link names not being propagated to this parser.
    """
    pw_calc = CalculationFactory('quantumespresso.pw')
    pw_spec = pw_calc.spec()
    assert 'structure' in pw_spec.inputs, list(pw_spec.inputs.keys())
    assert 'kpoints' in pw_spec.inputs, list(pw_spec.inputs.keys())
    assert 'output_parameters' in pw_spec.outputs, list(pw_spec.outputs.keys())
    assert 'output_structure' in pw_spec.outputs, list(pw_spec.outputs.keys())


def test_projwfc_default(
    fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser, parse_from_node,
    parser_fixture_path, projwfc_inputs, data_regression
):
    """Test ``ProjwfcParser`` on the results of a simple ``projwfc.x`` calculation."""
    entry_point_calc_job = 'quantumespresso.projwfc'
    entry_point_parser = 'quantumespresso.projwfc'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default', inputs=projwfc_inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parse_from_node(
        parser, node, store_provenance=False, retrieved_temp=os.path.join(parser_fixture_path, 'projwfc', 'default')
    )

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message

    for link_name in ['output_parameters', 'Dos', 'bands', 'projections', 'lowdin']:
        assert link_name in results, list(results.keys())

    data_regression.check({
        'Dos': results['Dos'].attributes,
        'bands': results['bands'].attributes,
        'projections':
        {k: v for k, v in results['projections'].attributes.items() if k not in ['reference_bandsdata_uuid']},
        'lowdin': {k: v for k, v in results['lowdin'].attributes.items() if k not in ['structure_uuid']}
    })
