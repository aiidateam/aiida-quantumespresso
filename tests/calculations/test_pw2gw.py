# -*- coding: utf-8 -*-
"""Tests for the `Pw2gwCalculation` class."""
from __future__ import absolute_import

from aiida import orm
from aiida.common import datastructures

from aiida_quantumespresso.utils.resources import get_default_options


def test_pw_default(
    aiida_profile, fixture_localhost, fixture_sandbox, fixture_code, generate_calc_job, generate_remote_data, tmpdir,
    file_regression
):
    """Test a default `Pw2gwCalculation`."""
    entry_point_name = 'quantumespresso.pw2gw'

    parameters = {
        'INPUTPP': {
            'qplda': False,
            'vxcdiag': False,
            'vkb': False,
            'Emin': 0.0,
            'Emax': 15.0,
            'DeltaE': 0.001,
        }
    }

    parent = generate_remote_data(
        fixture_localhost,
        str(tmpdir),
        'quantumespresso.pw',
    )

    inputs = {
        'code': fixture_code(entry_point_name),
        'parameters': orm.Dict(dict=parameters),
        'parent_folder': parent,
        'metadata': {
            'options': get_default_options()
        }
    }

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    retrieve_list = ['aiida.out', 'epsX.dat', 'epsY.dat', 'epsZ.dat', 'epsTOT.dat']

    assert isinstance(calc_info, datastructures.CalcInfo)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    with fixture_sandbox.open('aiida.in') as handle:
        input_written = handle.read()

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(['aiida.in'])
    file_regression.check(input_written, encoding='utf-8', extension='.in')
