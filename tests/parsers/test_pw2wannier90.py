# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `Pw2wannier90Parser`."""
from __future__ import absolute_import

import pytest
import os

from aiida import orm
from aiida.common import AttributeDict


@pytest.fixture
def pw2wannier90_inputs():
    """
    Minimal input for pw2wannier90 calculations.
    """
    basepath = os.path.dirname(os.path.abspath(__file__))
    nnkp_filepath = os.path.join(basepath, 'fixtures', 'pw2wannier90', 'inputs', 'aiida.nnkp')

    parameters = {
        'inputpp': {
            'write_amn': False,
            'write_mmn': False,
            'write_unk': False,
            'scdm_proj': True,
            'scdm_entanglement': 'isolated',
        }
    }

    settings = {'ADDITIONAL_RETRIEVE_LIST': ['*.amn', '*.mmn', '*.eig']}

    # Since we don't actually run pw2wannier.x, we only pretend to have the output folder
    # of a parent pw.x calculation. The nnkp file, instead, is real.
    inputs = {
        'parent_folder': orm.FolderData().store(),
        'nnkp_file': orm.SinglefileData(file=nnkp_filepath).store(),
        'parameters': orm.Dict(dict=parameters),
        'settings': orm.Dict(dict=settings),
    }

    return AttributeDict(inputs)


def test_pw2wannier90_default(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
                              pw2wannier90_inputs, data_regression):
    """
    Test a minimal `pw2wannier.x` calculation.
    The parser only checks for errors in aiida.out, so the reference contents of output_parameters
    will also be very minimal.
    """
    entry_point_calc_job = 'quantumespresso.pw2wannier90'
    entry_point_parser = 'quantumespresso.pw2wannier90'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default', pw2wannier90_inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    data_regression.check({
        'parameters': results['output_parameters'].get_dict()
    })
