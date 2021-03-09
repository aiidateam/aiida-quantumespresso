# -*- coding: utf-8 -*-
# pylint: disable=invalid-name,redefined-outer-name
"""Tests for the `PwParser`."""
import pytest

from aiida import orm
from aiida.common import AttributeDict


@pytest.fixture
def generate_inputs(generate_structure):
    """Return only those inputs that the parser will expect to be there."""

    def _generate_inputs(calculation_type='scf', parameters=None, settings=None, metadata=None):
        structure = generate_structure()
        parameters = {'CONTROL': {'calculation': calculation_type}, **(parameters or {})}
        kpoints = orm.KpointsData()
        kpoints.set_cell_from_structure(structure)
        kpoints.set_kpoints_mesh_from_density(0.15)

        return AttributeDict({
            'structure': generate_structure(),
            'kpoints': kpoints,
            'parameters': orm.Dict(dict=parameters),
            'settings': orm.Dict(dict=settings),
            'metadata': metadata or {}
        })

    return _generate_inputs


def test_pw_default(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression):
    """Test a `pw.x` calculation in `scf` mode.

    The output is created by running a dead simple SCF calculation for a silicon structure. This test should test the
    standard parsing of the stdout content and XML file stored in the standard results node.
    """
    name = 'default'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node), [log.message for log in orm.Log.objects.get_logs_for(node)]
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_trajectory' in results

    data_regression.check({
        'output_kpoints': results['output_kpoints'].attributes,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].attributes,
    })


def test_pw_default_no_xml(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test the parsing of an output directory without an XML file."""
    name = 'default_no_xml'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    # With default inputs, the parsing should fail
    inputs = generate_inputs(calculation_type='relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_failed, calcfunction.process_state
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_XML_MISSING.status

    # By setting the `without_xml` option the parsing should succeed and simply only use the stdout content
    inputs = generate_inputs(calculation_type='relax', metadata={'options': {'without_xml': True}})
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results

    data_regression.check({
        'output_parameters': results['output_parameters'].get_dict(),
        'output_structure': results['output_structure'].attributes,
        'output_trajectory': results['output_trajectory'].attributes,
    })


def test_pw_default_xml_200420(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test a `pw.x` calculation in `scf` mode that produced the XML output with schema of 200420.

    The output is created by running a dead simple SCF calculation for an aluminium structure. This test should test the
    standard parsing of the stdout content and XML file stored in the standard results node.
    """
    name = 'default_xml_200420'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node), [log.message for log in orm.Log.objects.get_logs_for(node)]
    assert 'output_band' in results
    assert 'output_parameters' in results
    assert 'output_trajectory' in results

    data_regression.check({
        'output_band': results['output_band'].attributes,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].attributes,
    })


def test_pw_default_xml_190304(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test a `pw.x` calculation in `scf` mode that produced the XML output with schema of 190304.

    The output is created by running a dead simple SCF calculation for a silicon structure. This test should test the
    standard parsing of the stdout content and XML file stored in the standard results node.
    """
    name = 'default_xml_190304'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node), [log.message for log in orm.Log.objects.get_logs_for(node)]
    assert 'output_band' in results
    assert 'output_parameters' in results
    assert 'output_trajectory' in results

    data_regression.check({
        'output_band': results['output_band'].attributes,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].attributes,
    })


def test_pw_default_xml_191206(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test a `pw.x` calculation in `scf` mode that produced the XML output with schema of 191206.

    The output is created by running a dead simple SCF calculation for an aluminium structure. This test should test the
    standard parsing of the stdout content and XML file stored in the standard results node.
    """
    name = 'default_xml_191206'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node), [log.message for log in orm.Log.objects.get_logs_for(node)]
    assert 'output_band' in results
    assert 'output_parameters' in results
    assert 'output_trajectory' in results

    data_regression.check({
        'output_band': results['output_band'].attributes,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].attributes,
    })


def test_pw_initialization_xml_new(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test a `pw.x` calculation with new XML that only runs the preamble, i.e. an initialization-only calculation."""
    name = 'initialization_xml_new'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(settings={'ONLY_INITIALIZATION': True})
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node), [log.message for log in orm.Log.objects.get_logs_for(node)]
    assert 'output_band' not in results
    assert 'output_kpoints' not in results
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].attributes,
    })


def test_pw_failed_computing_cholesky(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test the parsing of a calculation that failed during cholesky factorization.

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = 'failed_computing_cholesky'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_COMPUTING_CHOLESKY.status


def test_pw_failed_dexx_negative(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test the parsing of a calculation that failed due to negative dexx.

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = 'failed_dexx_negative'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_DEXX_IS_NEGATIVE.status


def test_pw_failed_missing(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test the parsing of a calculation that was interrupted before output files could even be written.

    In this particular interrupted test both the XML and the stdout are completely missing.

    This test simulates where a calculation fails to write output files entirely, probably due to grave crashes such
    as segmentation faults.
    """
    name = 'failed_missing'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_FILES.status
    assert orm.Log.objects.get_logs_for(node)


def test_pw_failed_interrupted(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
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

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_FILES.status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    data_regression.check(results['output_parameters'].get_dict())


def test_pw_failed_interrupted_stdout(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
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

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE.status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check(results['output_parameters'].get_dict())


def test_pw_failed_interrupted_xml(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
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

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_XML_PARSE.status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check(results['output_parameters'].get_dict())


def test_pw_npools_too_high(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test the parsing of a calculation that failed because some nodes have no k-points.

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = 'failed_npools_too_high'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_NPOOLS_TOO_HIGH.status


def test_tot_magnetization(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test the parsing of a calculation that specifies tot_magnetization.

    In this case there are two Fermi energies, see:

    https://lists.quantum-espresso.org/pipermail/users/2011-April/020089.html
    """
    name = 'tot_magnetization'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_status
    assert 'output_parameters' in results
    output_parameters = results['output_parameters'].get_dict()
    data_regression.check({
        'fermi_energy_up': output_parameters['fermi_energy_up'],
        'fermi_energy_down': output_parameters['fermi_energy_down']
    })


def test_pw_failed_out_of_walltime(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test the parsing of an scf calculation that ran nominally but was cut short because it ran out of walltime."""
    name = 'failed_out_of_walltime'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUT_OF_WALLTIME.status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].attributes,
    })


def test_pw_failed_out_of_walltime_interrupted(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test the parsing of an scf calculation that ran nominally but was cut short because it ran out of walltime.

    This differs from `test_pw_failed_out_of_walltime` in the sense that even though QE initiated the termination of the
    calculation due to the walltime being exceeded, before it could write all necessary files too disk, the scheduler
    killed the job because the walltime was exceeded.
    """
    name = 'failed_out_of_walltime_interrupted'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUT_OF_WALLTIME_INTERRUPTED.status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].attributes,
    })


def test_pw_failed_scf_not_converged(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test the parsing of an scf calculation that ran nominally but did not reach convergence."""
    name = 'failed_scf_not_converged'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED.status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].attributes,
    })


def test_pw_relax_success(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression):
    """Test a `relax` that successfully converges."""
    name = 'relax_success'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node), [log.message for log in orm.Log.objects.get_logs_for(node)]
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


def test_pw_relax_failed_electronic(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test a `relax` that failed to converge during electronic cycle before ionic convergence is reached."""
    name = 'relax_failed_electronic'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_ELECTRONIC_CONVERGENCE_NOT_REACHED.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results


def test_pw_relax_failed_not_converged_nstep(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs
):
    """Test a `relax` that failed to converge within the maximum number of steps."""
    name = 'relax_failed_not_converged_nstep'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_EXCEEDED_NSTEP.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results


def test_pw_vcrelax_success(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test a `vc-relax` that successfully converges and the final scf also converges."""
    name = 'vcrelax_success'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node), [log.message for log in orm.Log.objects.get_logs_for(node)]
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


def test_pw_vcrelax_success_fractional(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test a `vc-relax`, that successfully converges and the final scf also converges.

    In this case the input atomic positions were defined using 'crystal' (i.e. fractional) units.
    """
    name = 'vcrelax_success_fractional'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node), [log.message for log in orm.Log.objects.get_logs_for(node)]
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


def test_pw_vcrelax_success_rVV10(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test a `vc-relax` rVV10 run that successfully converges."""
    name = 'vcrelax_success_rVV10'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node), [log.message for log in orm.Log.objects.get_logs_for(node)]
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'energy_vdw': results['output_parameters']['energy_vdw'],
        'array|stress': results['output_trajectory'].attributes['array|stress'],
    })


def test_pw_vcrelax_success_external_pressure(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs
):
    """Test a `vc-relax` with external pressure that successfully converges and the final scf also converges."""
    name = 'vcrelax_success_external_pressure'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.objects.get_logs_for(node), [log.message for log in orm.Log.objects.get_logs_for(node)]
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results


def test_pw_vcrelax_success_atoms_shape(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test a relax calculation that fixes cell volume and only changes shape, that successfully converges.

    This is an example of a variable cell relaxation calculation where `CELL.cell_dofree` is set to anything other than
    the default `all` in which case the threshold on the stress/pressure should be ignored.
    """
    name = 'vcrelax_success_atoms_shape'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax', parameters={'CELL': {'cell_dofree': 'shape'}})
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results


def test_pw_vcrelax_failed_charge_wrong(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test a `vc-relax` that failed because the integrated charge is different from the expected one."""
    name = 'vcrelax_failed_charge_wrong'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_CHARGE_IS_WRONG.status

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == expected_exit_status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results


def test_pw_vcrelax_failed_symmetry_not_orthogonal(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs
):
    """Test a `vc-relax` that failed because original symmetries no longer map onto new structure."""
    name = 'vcrelax_failed_symmetry_not_orthogonal'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_SYMMETRY_NON_ORTHOGONAL_OPERATION.status

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == expected_exit_status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_parameters' in results


def test_pw_vcrelax_failed_bfgs_history(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test a `vc-relax` that failed to converge due to two consecutive failures of BFGS."""
    name = 'vcrelax_failed_bfgs_history'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results


def test_pw_vcrelax_failed_bfgs_history_already_converged(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs
):
    """Test a `vc-relax` that stops due to two consecutive failures of BFGS but is actually converged.

    Quantum ESPRESSO can stop with the BFGS history reset error even when all forces and stresses are already below the
    thresholds. This happens when the structure is already relaxed at the beginning of the calculation. For some reason
    QE will still enter the optimization loop and often BFGS will struggle and eventually fail and stop. In this case,
    the parser should just consider the structure as relaxed and return `0`.
    """
    name = 'vcrelax_failed_bfgs_history_already_converged'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished_ok, calcfunction.exit_status
    assert not orm.Log.objects.get_logs_for(node), [log.message for log in orm.Log.objects.get_logs_for(node)]
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results


def test_pw_vcrelax_failed_bfgs_history_final_scf(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs
):
    """Test a `vc-relax` that failed to converge due to two consecutive failures of BFGS and final SCF fails."""
    name = 'vcrelax_failed_bfgs_history_final_scf'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results


def test_pw_vcrelax_failed_electronic(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test a `vc-relax` that failed to converge during electronic cycle before ionic convergence is reached."""
    name = 'vcrelax_failed_electronic'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_ELECTRONIC_CONVERGENCE_NOT_REACHED.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results


def test_pw_vcrelax_failed_electronic_final_scf(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs
):
    """Test a `vc-relax` that failed to converge in electronic cycle in the final SCF after ionic convergence."""
    name = 'vcrelax_failed_electronic_final_scf'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_FINAL_SCF_FAILED.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results


def test_pw_vcrelax_failed_not_converged_final_scf(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs
):
    """Test a `vc-relax` that successfully converges in ionic cycle, but thresholds are exceeded in the SCF."""
    name = 'vcrelax_failed_not_converged_final_scf'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results


def test_pw_vcrelax_failed_not_converged_nstep(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs
):
    """Test a `vc-relax` that failed to converge within the maximum number of steps."""
    name = 'vcrelax_failed_not_converged_nstep'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_IONIC_CYCLE_EXCEEDED_NSTEP.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert orm.Log.objects.get_logs_for(node)
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
