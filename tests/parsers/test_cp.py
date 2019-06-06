# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `CpParser`."""
from __future__ import absolute_import

from aiida import orm
from aiida.common import AttributeDict
from ase.spacegroup import crystal
import pytest


@pytest.fixture
def generate_inputs():
    """Return only those inputs that the parser will expect to be there."""
    alat = 5.4
    ase = crystal(
        'Si',
        [(0, 0, 0)],
        spacegroup=227,
        cellpar=[alat, alat, alat, 90, 90, 90],
        primitive_cell=True,
    )
    structure = orm.StructureData(ase=ase)
    parameters = {
        'CONTROL': {
            'calculation': 'cp',
            'restart_mode': 'from_scratch',
            'wf_collect': False,
            'iprint': 1,
            'isave': 100,
            'dt': 3.0,
            'max_seconds': 25 * 60,
            'nstep': 10,
        },
        'SYSTEM': {
            'ecutwfc': 30.0,
            'ecutrho': 240.0,
            'nr1b': 24,
            'nr2b': 24,
            'nr3b': 24,
        },
        'ELECTRONS': {
            'electron_damping': 1.0e-1,
            'electron_dynamics': 'damp',
            'emass': 400.0,
            'emass_cutoff': 3.0,
        },
        'IONS': {
            'ion_dynamics': 'none'
        },
    }
    return AttributeDict({
        'structure': structure,
        'parameters': orm.Dict(dict=parameters),
    })


def test_cp_default(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
                    generate_inputs, data_regression):
    """Test a default `cp.x` calculation."""
    entry_point_calc_job = 'quantumespresso.cp'
    entry_point_parser = 'quantumespresso.cp'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default', generate_inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'parameters': results['output_parameters'].get_dict(),
        'trajectory': results['output_trajectory'].attributes
    })
