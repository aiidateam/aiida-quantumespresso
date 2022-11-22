# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the ``OpenGridParser``."""


def test_opengrid_default(fixture_localhost, generate_calc_job_node, generate_parser, data_regression, num_regression):
    """Test ``OpenGridParser`` on the results of a simple ``open_grid.x`` calculation."""
    entry_point_calc_job = 'quantumespresso.open_grid'
    entry_point_parser = 'quantumespresso.open_grid'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default')
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message

    for link_name in ['kpoints_mesh', 'kpoints']:
        assert link_name in results, list(results.keys())

    data_regression.check({
        'kpoints_mesh': results['kpoints_mesh'].attributes,
        'kpoints': results['kpoints'].attributes,
    })

    num_regression.check({
        'kpoints_array0': results['kpoints'].get_array('kpoints')[:, 0],
        'kpoints_array1': results['kpoints'].get_array('kpoints')[:, 1],
        'kpoints_array2': results['kpoints'].get_array('kpoints')[:, 2],
    })


def test_opengrid_fftgrid(fixture_localhost, generate_calc_job_node, generate_parser):
    """Test ``OpenGridParser`` on the parsing of 'incompatible FFT grid' error."""
    entry_point_calc_job = 'quantumespresso.open_grid'
    entry_point_parser = 'quantumespresso.open_grid'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'fftgrid')
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_INCOMPATIBLE_FFT_GRID.status
