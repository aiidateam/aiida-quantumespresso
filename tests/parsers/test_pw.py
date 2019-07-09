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
    structure.append_atom(position=(0., 0., 0.), symbols='Si', name="Si1")
    structure.append_atom(position=(a / 4., a / 4., a / 4.), symbols='Si', name="Si2")
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
    name = 'default'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_default)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_array' in results
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_trajectory' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive
    data_regression.check({
        'output_array': results['output_array'].attributes,
        'output_kpoints': results['output_kpoints'].attributes,
        'output_parameters': results['output_parameters'].get_dict()
    })


def test_pw_default_xml_new(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
        generate_inputs_default, data_regression):
    """Test a `pw.x` calculation in `scf` mode that produced the XML output in the new schema-based format.

    The output is created by running a dead simple SCF calculation for a silicon structure.
    This test should test the standard parsing of the stdout content and XML file stored in the standard results node.
    """
    name = 'default_xml_new'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_default)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_array' in results
    assert 'output_band' in results
    assert 'output_parameters' in results
    assert 'output_trajectory' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive
    data_regression.check({
        'output_array': results['output_array'].attributes,
        'output_band': results['output_band'].attributes,
        'output_parameters': results['output_parameters'].get_dict()
    })


def test_pw_failed_missing(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_default, data_regression):
    """Test the parsing of a calculation that was interrupted before output files could even be written.

    In this particular interrupted test both the XML and the stdout are completely missing.

    This test simulates where a calculation fails to write output files entirely, probably due to grave crashes such
    as segmentation faults.
    """
    name = 'failed_missing'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_default)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_FILES.status


def test_pw_failed_interrupted(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_default, data_regression):
    """Test the parsing of a calculation that was interrupted *after* convergence was achieved.

    In this particular interrupted test both the XML and the stdout are incomplete.

    This test simulates where an SCF calculation reaches convergence but the code is interrupted while writing the
    final output to disk. This can occur for a variety of reasons, for example the scheduler killing the job short
    due to out of walltime or out of memory errors.

    Only the output parameters are expected for the outputs since the `array` and `kpoints` are parsed from the XML.
    """
    name = 'failed_interrupted'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_default)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_FILES.status
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())


def test_pw_failed_interrupted_stdout(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_default, data_regression):
    """Test the parsing of a calculation that was interrupted *after* convergence was achieved.

    In this particular interrupted test only the stdout is incomplete and the XML is valid.

    This test simulates where an SCF calculation reaches convergence but the code is interrupted while writing the
    final output to disk. This can occur for a variety of reasons, for example the scheduler killing the job short
    due to out of walltime or out of memory errors.

    All three base outputs `array`, `kpoints` and `parameters` are expected as the first two are parsed from the XML
    which is in tact and the parameters are parsed from `stdout`, which, although interrupted, is mostly complete.
    """
    name = 'failed_interrupted_stdout'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_default)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE.status
    assert 'output_array' in results
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())


def test_pw_failed_interrupted_xml(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_default, data_regression):
    """Test the parsing of a calculation that was interrupted *after* convergence was achieved.

    In this particular interrupted test only the XML is incomplete and the stdout is valid.

    This test simulates where an SCF calculation reaches convergence but the code is interrupted while writing the
    final output to disk. This can occur for a variety of reasons, for example the scheduler killing the job short
    due to out of walltime or out of memory errors.

    Only the `kpoints` are not expected in the outputs, since it is parsed from the XML which is corrupted in this test.
    """
    name = 'failed_interrupted_xml'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_default)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_XML_PARSE.status
    assert 'output_array' in results
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())


def test_pw_failed_out_of_walltime(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_default, data_regression):
    """Test the parsing of an scf calculation that ran nominally but was cut short because it ran out of walltime."""
    name = 'failed_out_of_walltime'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_default)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUT_OF_WALLTIME.status
    assert 'output_array' in results
    assert 'output_parameters' in results
    data_regression.check({
        'output_array': results['output_array'].attributes,
        'output_parameters': results['output_parameters'].get_dict(),
    })


def test_pw_failed_scf_not_converged(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_default, data_regression):
    """Test the parsing of an scf calculation that ran nominally but did not reach convergence."""
    name = 'failed_scf_not_converged'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_default)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED.status
    assert 'output_array' in results
    assert 'output_parameters' in results
    data_regression.check({
        'output_array': results['output_array'].attributes,
        'output_parameters': results['output_parameters'].get_dict(),
    })


def test_pw_relax_success(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
        generate_inputs_relax, data_regression):
    """Test a `relax` that successfully converges."""
    name = 'relax_success'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_relax)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    assert 'output_array' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive
    data_regression.check({
        'output_kpoints': results['output_kpoints'].attributes,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_structure': results['output_structure'].attributes,
        'output_trajectory': results['output_trajectory'].attributes,
    })


def test_pw_relax_failed_electronic(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_relax, data_regression):
    """Test a `relax` that failed to converge during electronic cycle before ionic convergence is reached."""
    name = 'relax_failed_electronic'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_relax)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_ELECTRONIC_CONVERGENCE_NOT_REACHED.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    assert 'output_array' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive


def test_pw_relax_failed_not_converged_nstep(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_relax, data_regression):
    """Test a `relax` that failed to converge within the maximum number of steps."""
    name = 'relax_failed_not_converged_nstep'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_relax)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_EXCEEDED_NSTEP.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    assert 'output_array' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive


def test_pw_vcrelax_success(fixture_database, fixture_computer_localhost, generate_calc_job_node, generate_parser,
        generate_inputs_relax, data_regression):
    """Test a `vc-relax` that successfully converges and the final scf also converges."""
    name = 'vcrelax_success'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_relax)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    assert 'output_array' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive
    data_regression.check({
        'output_kpoints': results['output_kpoints'].attributes,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_structure': results['output_structure'].attributes,
        'output_trajectory': results['output_trajectory'].attributes,
    })


def test_pw_vcrelax_failed_bfgs_history(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_relax, data_regression):
    """Test a `vc-relax` that failed to converge due to two consecutive failures of BFGS."""
    name = 'vcrelax_failed_bfgs_history'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_relax)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    assert 'output_array' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive


def test_pw_vcrelax_failed_bfgs_history_final_scf(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_relax, data_regression):
    """Test a `vc-relax` that failed to converge due to two consecutive failures of BFGS and final SCF fails."""
    name = 'vcrelax_failed_bfgs_history_final_scf'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_relax)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    assert 'output_array' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive


def test_pw_vcrelax_failed_electronic(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_relax, data_regression):
    """Test a `vc-relax` that failed to converge during electronic cycle before ionic convergence is reached."""
    name = 'vcrelax_failed_electronic'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_relax)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_ELECTRONIC_CONVERGENCE_NOT_REACHED.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    assert 'output_array' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive


def test_pw_vcrelax_failed_electronic_final_scf(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_relax, data_regression):
    """Test a `vc-relax` that failed to converge in electronic cycle in the final SCF after ionic convergence."""
    name = 'vcrelax_failed_electronic_final_scf'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_relax)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_FINAL_SCF_FAILED.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    assert 'output_array' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive


def test_pw_vcrelax_failed_not_converged_final_scf(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_relax, data_regression):
    """Test a `vc-relax` that successfully converges in ionic cycle, but thresholds are exceeded in the SCF."""
    name = 'vcrelax_failed_not_converged_final_scf'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_relax)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    assert 'output_array' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive


def test_pw_vcrelax_failed_not_converged_nstep(fixture_database, fixture_computer_localhost, generate_calc_job_node,
        generate_parser, generate_inputs_relax, data_regression):
    """Test a `vc-relax` that failed to converge within the maximum number of steps."""
    name = 'vcrelax_failed_not_converged_nstep'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_computer_localhost, name, generate_inputs_relax)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_EXCEEDED_NSTEP.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    assert 'output_array' not in results  # The `ArrayData` and `TrajectoryData` outputs are mutually exclusive
