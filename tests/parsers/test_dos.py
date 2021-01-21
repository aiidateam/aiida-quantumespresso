# -*- coding: utf-8 -*-
"""Tests for the `DosParser`."""

from aiida import orm
from aiida.common import AttributeDict


def generate_inputs():
    """Return only those inputs that the parser will expect to be there."""
    return AttributeDict()


def test_dos_default(fixture_localhost, generate_calc_job_node, generate_parser, data_regression, num_regression):
    """Test `DosParser` on the results of a simple `dos.x` calculation."""
    entry_point_calc_job = 'quantumespresso.dos'
    entry_point_parser = 'quantumespresso.dos'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default', generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    out_params = results['output_parameters'].get_dict()
    out_dos_x = results['output_dos'].get_x()
    out_dos_y = results['output_dos'].get_y()
    dos_labels = [out_dos_x[0]] + [arr[0] for arr in out_dos_y]  # 4 strings
    dos_values = [out_dos_x[1]] + [arr[1] for arr in out_dos_y]  # 4 numpy arrays
    dos_units = [out_dos_x[2]] + [arr[2] for arr in out_dos_y]  # 4 strings

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_dos' in results
    data_regression.check({
        'parameters': out_params,
        'dos': {
            'labels': dos_labels,
            'units': dos_units,
        }
    })
    num_regression.check({f'dos_val_{i}': val for i, val in enumerate(dos_values)},
                         default_tolerance=dict(atol=0, rtol=1e-18))


def test_dos_failed_interrupted(fixture_localhost, generate_calc_job_node, generate_parser):
    """Test `DosParser` on the results of a `dos.x` calculation that was interrupted abruptly."""
    entry_point_calc_job = 'quantumespresso.dos'
    entry_point_parser = 'quantumespresso.dos'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'failed_interrupted', generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_failed
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE.status
