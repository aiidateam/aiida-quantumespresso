# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the `ProjwfcParser`."""

import pytest


@pytest.fixture
def generate_projwfc_node(generate_calc_job_node, fixture_localhost, tmpdir):
    """Fixture to constructure a ``projwfc.x`` calcjob node for a specified test."""

    def _generate_projwfc_node(test_name):
        """Generate a mock ``ProjwfcCalculation`` node for testing the parsing.

        :param test_name: The name of the test folder that contains the output files.
        """
        entry_point_calc_job = 'quantumespresso.projwfc'

        retrieve_temporary_list = ['data-file-schema.xml', '*.pdos*']
        attributes = {'retrieve_temporary_list': retrieve_temporary_list}

        node = generate_calc_job_node(
            entry_point_name=entry_point_calc_job,
            computer=fixture_localhost,
            test_name=test_name,
            attributes=attributes,
            retrieve_temporary=(tmpdir, retrieve_temporary_list)
        )
        return node

    return _generate_projwfc_node


@pytest.mark.parametrize('test_name', ('nonpolarised', 'noncollinear', 'spinorbit', 'numbered_kinds'))
def test_projwfc(generate_projwfc_node, generate_parser, data_regression, tmpdir, test_name):
    """Test ``ProjwfcParser`` on the results of a non-polarised ``projwfc.x`` calculation."""
    node = generate_projwfc_node(test_name)
    parser = generate_parser('quantumespresso.projwfc')
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message

    for link_name in ['output_parameters', 'Dos', 'bands', 'projections']:
        assert link_name in results, list(results.keys())

    data_regression.check({
        'Dos': results['Dos'].base.attributes.all,
        'bands': results['bands'].base.attributes.all,
        'projections':
        {k: v for k, v in results['projections'].base.attributes.all.items() if k not in ['reference_bandsdata_uuid']}
    })


def test_projwfc_spinpolarised(generate_projwfc_node, generate_parser, data_regression, tmpdir):
    """Test ``ProjwfcParser`` on the results of a spin-polarised ``projwfc.x`` calculation."""
    node = generate_projwfc_node('spinpolarised')
    parser = generate_parser('quantumespresso.projwfc')
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message

    for link_name in ['output_parameters', 'Dos', 'bands_up', 'bands_down', 'projections_up', 'projections_down']:
        assert link_name in results, list(results.keys())

    data_regression.check({
        'Dos':
        {array_name: results['Dos'].get_array(array_name).tolist() for array_name in results['Dos'].get_arraynames()},
        'bands_up': results['bands_up'].base.attributes.all,
        'bands_down': results['bands_down'].base.attributes.all,
        'projections_up': {
            k: v
            for k, v in results['projections_up'].base.attributes.all.items()
            if k not in ['reference_bandsdata_uuid']
        },
        'projections_down': {
            k: v
            for k, v in results['projections_down'].base.attributes.all.items()
            if k not in ['reference_bandsdata_uuid']
        }
    })


def test_projwfc_no_retrieved_temporary(generate_calc_job_node, fixture_localhost, generate_parser):
    """Test ``ProjwfcParser`` fails when the retrieved temporary folder is missing."""
    node = generate_calc_job_node(
        entry_point_name='quantumespresso.projwfc',
        computer=fixture_localhost,
        test_name='xml_missing',
    )
    parser = generate_parser('quantumespresso.projwfc')
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_failed, calcfunction.process_state
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_NO_RETRIEVED_TEMPORARY_FOLDER.status


@pytest.mark.parametrize('test_name, exit_status', (
    ['xml_missing', 303],
    ['xml_parse', 321],
    ['xml_format', 322],
))
def test_projwfc_xml_failures(generate_projwfc_node, generate_parser, tmpdir, test_name, exit_status):
    """Test ``ProjwfcParser`` fails when the XML is missing."""
    node = generate_projwfc_node(test_name)
    parser = generate_parser('quantumespresso.projwfc')
    _, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

    assert calcfunction.is_failed, calcfunction.process_state
    assert calcfunction.exit_status == exit_status
