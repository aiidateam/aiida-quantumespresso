# -*- coding: utf-8 -*-
"""Tests for the `XspectraParser`."""

from aiida import orm
from aiida.common import AttributeDict


def generate_inputs():
    """Return only those inputs that the parser will expect to be there."""
    return AttributeDict()


def test_xspectra_default(fixture_localhost, generate_calc_job_node, generate_parser, data_regression, num_regression):
    """Test `XspectraParser` on the results of a simple `xspectra.x` calculation."""
    entry_point_calc_job = 'quantumespresso.xspectra'
    entry_point_parser = 'quantumespresso.xspectra'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default', generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    out_params = results['output_parameters'].get_dict()
    out_spectra_x = results['spectra'].get_x()
    out_spectra_y = results['spectra'].get_y()
    spectra_labels = [out_spectra_x[0]] + [arr[0] for arr in out_spectra_y]  # 2 strings
    spectra_values = [out_spectra_x[1]] + [arr[1] for arr in out_spectra_y]  # 2 numpy arrays
    spectra_units = [out_spectra_x[2]] + [arr[2] for arr in out_spectra_y]  # 2 strings

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'spectra' in results
    data_regression.check({
        'parameters': out_params,
        'xspectra': {
            'labels': spectra_labels,
            'units': spectra_units,
        }
    })
    num_regression.check({f'spectra_val_{i}': val for i, val in enumerate(spectra_values)},
                         default_tolerance=dict(atol=0, rtol=1e-10))


def test_xspectra_spin(fixture_localhost, generate_calc_job_node, generate_parser, data_regression, num_regression):
    """Test `XspectraParser` on the results of a spin-polarised `xspectra.x` calculation."""
    entry_point_calc_job = 'quantumespresso.xspectra'
    entry_point_parser = 'quantumespresso.xspectra'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'spin', generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    out_params = results['output_parameters'].get_dict()
    out_spectra_x = results['spectra'].get_x()
    out_spectra_y = results['spectra'].get_y()
    spectra_labels = [out_spectra_x[0]] + [arr[0] for arr in out_spectra_y]  # 4 strings
    spectra_values = [out_spectra_x[1]] + [arr[1] for arr in out_spectra_y]  # 4 numpy arrays
    spectra_units = [out_spectra_x[2]] + [arr[2] for arr in out_spectra_y]  # 4 strings

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'spectra' in results
    data_regression.check({
        'parameters': out_params,
        'xspectra': {
            'labels': spectra_labels,
            'units': spectra_units,
        }
    })
    num_regression.check({f'spectra_val_{i}': val for i, val in enumerate(spectra_values)},
                         default_tolerance=dict(atol=0, rtol=1e-10))


def test_xspectra_failed_interrupted(fixture_localhost, generate_calc_job_node, generate_parser):
    """Test `XspectraParser` on the results of an `xspectra.x` calculation that was interrupted abruptly."""
    entry_point_calc_job = 'quantumespresso.xspectra'
    entry_point_parser = 'quantumespresso.xspectra'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'failed_interrupted', generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_failed
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE.status


def test_xspectra_failed_xiabs_wrong(fixture_localhost, generate_calc_job_node, generate_parser):
    """Test `XspectraParser` to correctly report error code 313 (xiabs set to wrong element)."""
    entry_point_calc_job = 'quantumespresso.xspectra'
    entry_point_parser = 'quantumespresso.xspectra'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'failed_xiabs_wrong', generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_failed
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_ABSORBING_SPECIES_WRONG.status


def test_xspectra_failed_xiabs_zero(fixture_localhost, generate_calc_job_node, generate_parser):
    """Test `XspectraParser` to correctly report error code 314 (xiabs<1 or xiabs>ntyp)."""
    entry_point_calc_job = 'quantumespresso.xspectra'
    entry_point_parser = 'quantumespresso.xspectra'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'failed_xiabs_zero', generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_failed
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_ABSORBING_SPECIES_ZERO.status
