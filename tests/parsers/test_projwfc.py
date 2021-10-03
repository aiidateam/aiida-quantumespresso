# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the `ProjwfcParser`."""


def test_projwfc_nonpolarised(fixture_localhost, generate_calc_job_node, generate_parser, data_regression, tmpdir):
    """Test ``ProjwfcParser`` on the results of a non-polarised ``projwfc.x`` calculation."""
    entry_point_calc_job = 'quantumespresso.projwfc'
    entry_point_parser = 'quantumespresso.projwfc'

    retrieve_temporary_list = ['data-file-schema.xml']
    attributes = {'retrieve_temporary_list': retrieve_temporary_list}

    node = generate_calc_job_node(
        entry_point_name=entry_point_calc_job,
        computer=fixture_localhost,
        test_name='nonpolarised',
        attributes=attributes,
        retrieve_temporary=(tmpdir, retrieve_temporary_list)
    )
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

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


def test_projwfc_spinpolarised(fixture_localhost, generate_calc_job_node, generate_parser, data_regression, tmpdir):
    """Test ``ProjwfcParser`` on the results of a spin-polarised ``projwfc.x`` calculation."""
    entry_point_calc_job = 'quantumespresso.projwfc'
    entry_point_parser = 'quantumespresso.projwfc'

    retrieve_temporary_list = ['data-file-schema.xml']
    attributes = {'retrieve_temporary_list': retrieve_temporary_list}

    node = generate_calc_job_node(
        entry_point_name=entry_point_calc_job,
        computer=fixture_localhost,
        test_name='spinpolarised',
        attributes=attributes,
        retrieve_temporary=(tmpdir, retrieve_temporary_list)
    )
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

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


def test_projwfc_noncollinear(fixture_localhost, generate_calc_job_node, generate_parser, data_regression, tmpdir):
    """Test ``ProjwfcParser`` on the results of a noncollinear ``projwfc.x`` calculation."""
    entry_point_calc_job = 'quantumespresso.projwfc'
    entry_point_parser = 'quantumespresso.projwfc'

    retrieve_temporary_list = ['data-file-schema.xml']
    attributes = {'retrieve_temporary_list': retrieve_temporary_list}

    node = generate_calc_job_node(
        entry_point_name=entry_point_calc_job,
        computer=fixture_localhost,
        test_name='noncollinear',
        attributes=attributes,
        retrieve_temporary=(tmpdir, retrieve_temporary_list)
    )
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

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


def test_projwfc_spinorbit(fixture_localhost, generate_calc_job_node, generate_parser, data_regression, tmpdir):
    """Test ``ProjwfcParser`` on the results of a spinorbit ``projwfc.x`` calculation."""
    entry_point_calc_job = 'quantumespresso.projwfc'
    entry_point_parser = 'quantumespresso.projwfc'

    retrieve_temporary_list = ['data-file-schema.xml']
    attributes = {'retrieve_temporary_list': retrieve_temporary_list}

    node = generate_calc_job_node(
        entry_point_name=entry_point_calc_job,
        computer=fixture_localhost,
        test_name='spinorbit',
        attributes=attributes,
        retrieve_temporary=(tmpdir, retrieve_temporary_list)
    )
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

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
