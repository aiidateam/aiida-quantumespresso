# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `PwParser`."""
from __future__ import absolute_import

import pytest


def test_pw_default(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test a default `pw.x` calculation.

    The output is created by running a dead simple SCF calculation for a silicon structure.
    This test should test the standard parsing of the stdout content and XML file stored in the standard results node.
    """
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    input_nodes = generate_inputs(entry_point_calc_job)

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default', input_nodes)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results

    output_parameters = results['output_parameters']

    assert output_parameters['beta_real_space'] is False
    assert output_parameters['charge_density'] == './charge-density.dat'
    assert output_parameters['constraint_mag'] == 0
    assert output_parameters['creator_name'] == 'pwscf'
    assert output_parameters['creator_version'] == '6.1'
    assert output_parameters['dft_exchange_correlation'] == 'PBE'
    assert output_parameters['do_not_use_time_reversal'] is False
    assert output_parameters['energy'] == pytest.approx(-308.191875414092)
    assert output_parameters['energy_accuracy'] == pytest.approx(7.347073531662e-06)
    assert output_parameters['energy_accuracy_units'] == 'eV'
    assert output_parameters['energy_ewald'] == pytest.approx(-228.561248748643)
    assert output_parameters['energy_ewald_units'] == 'eV'
    assert output_parameters['energy_hartree'] == pytest.approx(17.2680757695669)
    assert output_parameters['energy_hartree_units'] == 'eV'
    assert output_parameters['energy_one_electron'] == pytest.approx(71.7330877993463)
    assert output_parameters['energy_one_electron_units'] == 'eV'
    assert output_parameters['energy_threshold'] == pytest.approx(3.84e-06)
    assert output_parameters['energy_units'] == 'eV'
    assert output_parameters['energy_xc'] == pytest.approx(-168.631790234362)
    assert output_parameters['energy_xc_units'] == 'eV'
    assert output_parameters['fermi_energy'] == pytest.approx(6.5091589497145)
    assert output_parameters['fermi_energy_units'] == 'eV'
    assert output_parameters['fft_grid'] == [36, 36, 36]
    assert output_parameters['format_name'] == 'qexml'
    assert output_parameters['format_version'] == '1.4.0'
    assert output_parameters['has_dipole_correction'] is False
    assert output_parameters['has_electric_field'] is False
    assert output_parameters['init_wall_time_seconds'] == 0.9
    assert output_parameters['inversion_symmetry'] is True
    assert output_parameters['lda_plus_u_calculation'] is False
    assert output_parameters['lkpoint_dir'] is True
    assert output_parameters['lsda'] is False
    assert output_parameters['magnetization_angle1'] == [0.0]
    assert output_parameters['magnetization_angle2'] == [0.0]
    assert output_parameters['monkhorst_pack_grid'] == [2, 2, 2]
    assert output_parameters['monkhorst_pack_offset'] == [0, 0, 0]
    assert output_parameters['no_time_rev_operations'] is False
    assert output_parameters['non_colinear_calculation'] is False
    assert output_parameters['number_of_atomic_wfc'] == 8
    assert output_parameters['number_of_atoms'] == 2
    assert output_parameters['number_of_bands'] == 4
    assert output_parameters['number_of_bravais_symmetries'] == 48
    assert output_parameters['number_of_electrons'] == 8.0
    assert output_parameters['number_of_k_points'] == 3
    assert output_parameters['number_of_species'] == 1
    assert output_parameters['number_of_spin_components'] == 1
    assert output_parameters['number_of_symmetries'] == 48
    assert output_parameters['parser_info'] == 'aiida-quantumespresso parser pw.x v3.0.0a1'
    assert output_parameters['parser_version'] == '3.0.0a1'
    assert output_parameters['parser_warnings'] == []
    assert output_parameters['pp_check_flag'] is True
    assert output_parameters['q_real_space'] is False
    assert output_parameters['rho_cutoff'] == pytest.approx(3265.366014072)
    assert output_parameters['rho_cutoff_units'] == 'eV'
    assert output_parameters['scf_iterations'] == 5
    assert output_parameters['smooth_fft_grid'] == [25, 25, 25]
    assert output_parameters['spin_orbit_calculation'] is False
    assert output_parameters['spin_orbit_domag'] is False
    assert output_parameters['starting_magnetization'] == [0.0]
    assert output_parameters['symmetries_units'] == 'crystal'
    assert output_parameters['time_reversal_flag'] is True
    assert output_parameters['total_number_of_scf_iterations'] == 5
    assert output_parameters['volume'] == pytest.approx(40.02575175)
    assert output_parameters['wall_time'] == '         1.86s '
    assert output_parameters['wall_time_seconds'] == 1.86
    assert output_parameters['warnings'] == []
    assert output_parameters['wfc_cutoff'] == pytest.approx(408.170751759)
    assert output_parameters['wfc_cutoff_units'] == 'eV'
    assert output_parameters['symmetries'] == [
        {'symmetry_number': 0, 't_rev': '0'},
        {'symmetry_number': 1, 't_rev': '0'},
        {'symmetry_number': 2, 't_rev': '0'},
        {'symmetry_number': 3, 't_rev': '0'},
        {'symmetry_number': 4, 't_rev': '0'},
        {'symmetry_number': 5, 't_rev': '0'},
        {'symmetry_number': 6, 't_rev': '0'},
        {'symmetry_number': 7, 't_rev': '0'},
        {'symmetry_number': 8, 't_rev': '0'},
        {'symmetry_number': 9, 't_rev': '0'},
        {'symmetry_number': 10, 't_rev': '0'},
        {'symmetry_number': 11, 't_rev': '0'},
        {'symmetry_number': 12, 't_rev': '0'},
        {'symmetry_number': 13, 't_rev': '0'},
        {'symmetry_number': 14, 't_rev': '0'},
        {'symmetry_number': 15, 't_rev': '0'},
        {'symmetry_number': 16, 't_rev': '0'},
        {'symmetry_number': 17, 't_rev': '0'},
        {'symmetry_number': 18, 't_rev': '0'},
        {'symmetry_number': 19, 't_rev': '0'},
        {'symmetry_number': 20, 't_rev': '0'},
        {'symmetry_number': 21, 't_rev': '0'},
        {'symmetry_number': 22, 't_rev': '0'},
        {'symmetry_number': 23, 't_rev': '0'},
        {'symmetry_number': 32, 't_rev': '0'},
        {'symmetry_number': 33, 't_rev': '0'},
        {'symmetry_number': 34, 't_rev': '0'},
        {'symmetry_number': 35, 't_rev': '0'},
        {'symmetry_number': 36, 't_rev': '0'},
        {'symmetry_number': 37, 't_rev': '0'},
        {'symmetry_number': 38, 't_rev': '0'},
        {'symmetry_number': 39, 't_rev': '0'},
        {'symmetry_number': 40, 't_rev': '0'},
        {'symmetry_number': 41, 't_rev': '0'},
        {'symmetry_number': 42, 't_rev': '0'},
        {'symmetry_number': 43, 't_rev': '0'},
        {'symmetry_number': 44, 't_rev': '0'},
        {'symmetry_number': 45, 't_rev': '0'},
        {'symmetry_number': 46, 't_rev': '0'},
        {'symmetry_number': 47, 't_rev': '0'},
        {'symmetry_number': 48, 't_rev': '0'},
        {'symmetry_number': 49, 't_rev': '0'},
        {'symmetry_number': 50, 't_rev': '0'},
        {'symmetry_number': 51, 't_rev': '0'},
        {'symmetry_number': 52, 't_rev': '0'},
        {'symmetry_number': 53, 't_rev': '0'},
        {'symmetry_number': 54, 't_rev': '0'},
        {'symmetry_number': 55, 't_rev': '0'}
    ]
