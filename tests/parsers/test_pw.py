# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `PwParser`."""
from __future__ import absolute_import

import pytest
from aiida import orm
from aiida.common import AttributeDict


@pytest.fixture
def generate_inputs_default():
    """Return only those inputs that the parser will expect to be there."""
    structure = orm.StructureData()
    parameters = {
        'CONTROL': {
            'calculation': 'scf'
        },
        'SYSTEM': {
            'ecutrho': 240.0,
            'ecutwfc': 30.0
        }
    }
    kpoints = orm.KpointsData()
    kpoints.set_cell_from_structure(structure)
    kpoints.set_kpoints_mesh_from_density(0.15)

    return AttributeDict({
        'structure': structure,
        'kpoints': kpoints,
        'parameters': orm.Dict(dict=parameters),
        'settings': orm.Dict()
    })


@pytest.fixture
def generate_inputs_relax():
    """Return only those inputs that the parser will expect to be there.

    This needs a separate input generation function from the default one, because the parser depends on certain values
    in the input parameters to determine what kind of calculation it was. For example, it will check the card
    `CONTROL.calculation` to determine whether the `TrajectoryData` should be attached. If we would not set it to
    `relax`, the parser would not parse that output node and the test would fail. Until we can make the raw output
    parser independent of the input parameters, this will have to remain a separate test inputs generator.
    """
    a = 5.43
    structure = orm.StructureData(cell=[[a / 2., a / 2., 0], [a / 2., 0, a / 2.], [0, a / 2., a / 2.]])
    structure.append_atom(position=(0., 0., 0.), symbols='Si')
    structure.append_atom(position=(a / 4., a / 4., a / 4.), symbols='Si')
    structure.store()

    parameters = {
        'CONTROL': {
            'calculation': 'relax'
        },
        'SYSTEM': {
            'ecutrho': 240.0,
            'ecutwfc': 30.0
        }
    }
    kpoints = orm.KpointsData()
    kpoints.set_cell_from_structure(structure)
    kpoints.set_kpoints_mesh_from_density(0.15)

    return AttributeDict({
        'structure': structure,
        'kpoints': kpoints,
        'parameters': orm.Dict(dict=parameters),
        'settings': orm.Dict()
    })


def test_pw_default(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
                    generate_inputs_default, data_regression):
    """Test a `pw.x` calculation in `scf` mode.

    The output is created by running a dead simple SCF calculation for a silicon structure.
    This test should test the standard parsing of the stdout content and XML file stored in the standard results node.
    """
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default', generate_inputs_default)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())


def test_pw_default_xml_new(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
                    generate_inputs_default, data_regression):
    """Test a `pw.x` calculation in `scf` mode that produced the XML output in the new schema-based format.

    The output is created by running a dead simple SCF calculation for a silicon structure.
    This test should test the standard parsing of the stdout content and XML file stored in the standard results node.
    """
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'default_xml_new', generate_inputs_default)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())


def test_pw_relax(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
                  generate_inputs_relax, data_regression):
    """Test a `pw.x` calculation in `relax` mode.

    The output is created by running a dead 'relax' calculation for a silicon structure.
    Note that `output_band` will not be there, because the bands are only parsed if the `retrieved_temporary_folder` is
    passed containing the necessary files. Since we don't pass those in this test, the output node will not be created.
    """
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, 'relax', generate_inputs_relax)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_kpoints': results['output_kpoints'].attributes,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_structure': results['output_structure'].attributes,
        'output_trajectory': results['output_trajectory'].attributes,
    })
