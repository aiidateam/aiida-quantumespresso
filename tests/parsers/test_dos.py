# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `DosParser`."""
from __future__ import absolute_import
import pytest

from aiida import orm
from aiida.common import AttributeDict


@pytest.fixture
def dos_inputs():
    # The DosParser doesn't really need to access the parent folder, but we'd like to make the inputs as realistic
    # as possible, so we create an empty FolderData and attach it as an input to the current CalcJobNode.
    inputs = {
        'parent_folder': orm.FolderData().store(),
    }

    return AttributeDict(inputs)


def test_dos_default(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
                     dos_inputs, data_regression, num_regression):
    """
    Test `DosParser` on the results of a simple `dos.x` calculation.
    """
    entry_point_calc_job = 'quantumespresso.dos'
    entry_point_parser = 'quantumespresso.dos'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default', dos_inputs)
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
    num_regression.check({
        'dos_val_{}'.format(i): val for i, val in enumerate(dos_values)
    }, default_tolerance=dict(atol=0, rtol=1e-18))
