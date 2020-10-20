# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the `ProjwfcParser`."""
import pytest

from aiida import orm
from aiida.common import AttributeDict, LinkType


@pytest.fixture
def generate_inputs(generate_calc_job_node, fixture_localhost, generate_structure, generate_kpoints_mesh):
    """Create the required inputs for the ``ProjwfcCalculation``."""

    def _generate_inputs(output_parameters=None):
        entry_point_name = 'quantumespresso.pw'

        parent_inputs = {
            'structure': generate_structure(),
            'kpoints': generate_kpoints_mesh(4),
        }
        parent = generate_calc_job_node(entry_point_name, fixture_localhost, 'default', inputs=parent_inputs)
        output_parameters = output_parameters or orm.Dict(dict={'number_of_spin_components': 1})
        output_parameters.add_incoming(parent, link_type=LinkType.CREATE, link_label='output_parameters')
        output_parameters.store()

        inputs = {
            'parent_folder': parent.outputs.remote_folder,
        }

        return AttributeDict(inputs)

    return _generate_inputs


def test_projwfc_default(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression):
    """Test ``ProjwfcParser`` on the results of a simple ``projwfc.x`` calculation."""
    entry_point_calc_job = 'quantumespresso.projwfc'
    entry_point_parser = 'quantumespresso.projwfc'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default', generate_inputs())
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


def test_projwfc_spin(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression):
    """Test ``ProjwfcParser`` on the results of a LSDA ``projwfc.x`` calculation."""
    entry_point_calc_job = 'quantumespresso.projwfc'
    entry_point_parser = 'quantumespresso.projwfc'

    output_parameters = orm.Dict(dict={'number_of_spin_components': 2})
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'spin', generate_inputs(output_parameters))
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message

    for link_name in ['output_parameters', 'Dos', 'bands_up', 'bands_down', 'projections_up', 'projections_down']:
        assert link_name in results, list(results.keys())

    data_regression.check({
        'Dos': results['Dos'].attributes,
        'bands_up': results['bands_up'].attributes,
        'bands_down': results['bands_down'].attributes,
        'projections_up':
        {k: v for k, v in results['projections_up'].attributes.items() if k not in ['reference_bandsdata_uuid']},
        'projections_down':
        {k: v for k, v in results['projections_down'].attributes.items() if k not in ['reference_bandsdata_uuid']}
    })


def test_projwfc_noncollinear(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test ``ProjwfcParser`` on the results of a noncollinear ``projwfc.x`` calculation."""
    entry_point_calc_job = 'quantumespresso.projwfc'
    entry_point_parser = 'quantumespresso.projwfc'

    output_parameters = orm.Dict(dict={'number_of_spin_components': 4, 'non_colinear_calculation': True})
    node = generate_calc_job_node(
        entry_point_calc_job, fixture_localhost, 'noncollinear', generate_inputs(output_parameters)
    )
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


def test_projwfc_spinorbit(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test ``ProjwfcParser`` on the results of a spinorbit ``projwfc.x`` calculation."""
    entry_point_calc_job = 'quantumespresso.projwfc'
    entry_point_parser = 'quantumespresso.projwfc'

    output_parameters = orm.Dict(
        dict={
            'number_of_spin_components': 4,
            'non_colinear_calculation': True,
            'spin_orbit_calculation': True
        }
    )
    node = generate_calc_job_node(
        entry_point_calc_job, fixture_localhost, 'spinorbit', generate_inputs(output_parameters)
    )
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
