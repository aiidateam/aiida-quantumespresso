# -*- coding: utf-8 -*-
# pylint: disable=invalid-name,redefined-outer-name, too-many-lines
"""Tests for the `PwParser`."""
from aiida import orm
from aiida.common import AttributeDict
import pytest

from aiida_quantumespresso.calculations.pw import PwCalculation


@pytest.fixture
def generate_inputs(generate_structure):
    """Return only those inputs that the parser will expect to be there."""

    def _generate_inputs(calculation_type='scf', parameters=None, settings=None, metadata=None):
        structure = generate_structure()
        parameters = {'CONTROL': {'calculation': calculation_type}, **(parameters or {})}
        kpoints = orm.KpointsData()
        kpoints.set_cell_from_structure(structure)
        kpoints.set_kpoints_mesh_from_density(0.15)

        inputs = {
            'structure': generate_structure(),
            'kpoints': kpoints,
            'parameters': orm.Dict(parameters),
            'metadata': metadata or {}
        }
        if settings:
            inputs['settings'] = orm.Dict(settings)

        return AttributeDict(inputs)

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
    assert not orm.Log.collection.get_logs_for(node), [log.message for log in orm.Log.collection.get_logs_for(node)]
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_trajectory' in results

    data_regression.check({
        'output_kpoints': results['output_kpoints'].base.attributes.all,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].base.attributes.all,
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
        'output_structure': results['output_structure'].base.attributes.all,
        'output_trajectory': results['output_trajectory'].base.attributes.all,
    })


@pytest.mark.parametrize(
    'xml_format', [
        '190304',
        '191206',
        '200420',
        '210716',
        '211101',
        '220603',
        '230310',
        '240411',
        '241015',
    ]
)
def test_pw_default_xml(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression, xml_format
):
    """Test a `pw.x` calculation in `scf` mode that produced the XML output with the supported schemas.

    The output is created by running a dead simple SCF calculation for an aluminium structure. This test should test the
    standard parsing of the stdout content and XML file stored in the standard results node.
    """
    name = f'default_xml_{xml_format}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert not orm.Log.collection.get_logs_for(node), [log.message for log in orm.Log.collection.get_logs_for(node)]
    assert 'output_band' in results
    assert 'output_parameters' in results
    assert 'output_trajectory' in results

    data_regression.check({
        'output_band': results['output_band'].base.attributes.all,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].base.attributes.all,
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
    assert not orm.Log.collection.get_logs_for(node), [log.message for log in orm.Log.collection.get_logs_for(node)]
    assert 'output_band' not in results
    assert 'output_kpoints' not in results
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].base.attributes.all,
    })


def test_pw_failed_base_exception(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, monkeypatch
):
    """Test the parsing of a calculation where an unknown exception is raised.

    This should return ``ERROR_UNEXPECTED_PARSER_EXCEPTION`` formatted with exception title.
    """
    from aiida_quantumespresso.parsers.parse_raw import pw

    exception = 'the parser encountered an error.'

    def parse_xml(*_, **__):
        return {}, {}

    def parse_stdout(*_, **__):
        raise RuntimeError(exception)

    name = 'failed_base_exception'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)

    monkeypatch.setattr(pw, 'parse_stdout', parse_stdout)
    monkeypatch.setattr(parser, 'parse_xml', parse_xml)

    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION.status
    assert exception in calcfunction.exit_message


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_pw_failed_computing_cholesky(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename
):
    """Test the parsing of a calculation that failed during cholesky factorization.

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = f'failed_computing_cholesky{filename}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_COMPUTING_CHOLESKY.status


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_failed_too_many_bands_not_converged(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename
):
    """Test the parsing of a calculation that failed during cholesky factorization.

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = f'failed_too_many_bands_not_converged{filename}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    desired_exit_status = node.process_class.exit_codes.ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED.status
    assert calcfunction.exit_status == desired_exit_status


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_pw_failed_dexx_negative(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename):
    """Test the parsing of a calculation that failed due to negative dexx.

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = f'failed_dexx_negative{filename}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_DEXX_IS_NEGATIVE.status


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_failed_s_matrix_not_positive_definite(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename
):
    """Test the parsing of a calculation for which the overlap S matrix was not positive definitive (Davidson).

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = f'failed_s_matrix_not_positive_definite{filename}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE.status


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_failed_zhegvd(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename):
    """Test the parsing of a calculation for which the ``zhegvd`` failed (PPCG).

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = f'failed_zhegvd_failed{filename}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_ZHEGVD_FAILED.status


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_failed_qr(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename):
    """Test the parsing of a calculation for which the ``[Q, R] = qr(X, 0)``failed (PPCG).

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = f'failed_qr_failed{filename}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_QR_FAILED.status


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_failed_eigenvectors_convergence(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename
):
    """Test the parsing of a calculation that failed to converge the eigenvectors.

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = f'failed_eigenvectors_convergence{filename}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_EIGENVECTOR_CONVERGENCE.status


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_failed_broyden_factorization(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename
):
    """Test the parsing of a calculation that failed the factorization in the Broyden routine.

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = f'failed_broyden_factorization{filename}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_BROYDEN_FACTORIZATION.status


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_failed_g_par(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename):
    """Test the parsing of a calculation that failed to find unique G vector (finite electric field routine).

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = f'failed_g_par{filename}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_G_PAR.status


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
    assert orm.Log.collection.get_logs_for(node)


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
    assert orm.Log.collection.get_logs_for(node)
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
    assert orm.Log.collection.get_logs_for(node)
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
    assert orm.Log.collection.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check(results['output_parameters'].get_dict())


@pytest.mark.parametrize(
    'test_case, expected_exit_code', (
        ('default', None),
        ('failed_interrupted', PwCalculation.exit_codes.ERROR_SCHEDULER_OUT_OF_WALLTIME),
    )
)
def test_pw_failed_interrupted_scheduler(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, test_case, expected_exit_code
):
    """Test that an exit code set by the scheduler is not overridden unless a more specific error is parsed.

    The test is run twice, once for the ``default`` test case and once for ``failed_interrupted``, which correspond to a
    successful run and a run that got interrupted (usually due to scheduler killing the job). Before calling the parser
    the ``ERROR_SCHEDULER_OUT_OF_WALLTIME`` is set on the node. In the case of the ``default`` test case, this should be
    ignored and the parser should return ``ExitCode(0)``. For the interrupted case, the exit code of the scheduler
    should not be overridden so the parser should return ``ERROR_SCHEDULER_OUT_OF_WALLTIME``.
    """
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    # Generate the node and set an exit status as if it would have been set by the scheduler parser
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, test_case, generate_inputs())
    node.set_exit_status(PwCalculation.exit_codes.ERROR_SCHEDULER_OUT_OF_WALLTIME.status)

    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    if expected_exit_code is None:
        assert calcfunction.is_finished_ok, (calcfunction.exit_status, calcfunction.exception)
    else:
        assert not calcfunction.is_finished_ok, calcfunction.exception
        assert calcfunction.exit_status == expected_exit_code.status


def test_pw_failed_interrupted_relax(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test that a `relax` calculation that was interrupted by the scheduler is nevertheless (partially) parsed.

    In this case, since a partial trajectory could be parsed from which one can restart, the scheduler code is replaced
    with the more specific ``ERROR_IONIC_INTERRUPTED_PARTIAL_TRAJECTORY``.
    """
    name = 'failed_interrupted_relax'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    node.set_exit_status(PwCalculation.exit_codes.ERROR_SCHEDULER_OUT_OF_WALLTIME.status)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_failed
    assert calcfunction.exit_status == PwCalculation.exit_codes.ERROR_IONIC_INTERRUPTED_PARTIAL_TRAJECTORY.status
    assert orm.Log.collection.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_pw_npools_too_high_error(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename
):
    """Test the parsing of a calculation that failed because some nodes have no k-points.

    In this test the stdout is incomplete, and the XML is missing completely. The stdout contains
    the relevant error message.
    """
    name = f'failed_npools_too_high{filename}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_NPOOLS_TOO_HIGH.status


def test_pw_npools_too_high_not_error(fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs):
    """Test the parsing of a success calculation that report 'some nodes have no k-points'.

    For new QE version (test on v6.8) 'some nodes have no k-points' is not raised as an error and stop the
    calculation. The output is different but still contain the same content, so instead of check content in
    line use regex match.
    """
    name = 'finished_npools_too_high'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs())
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results


@pytest.mark.parametrize('calculation', ('relax', 'vc-relax'))
@pytest.mark.parametrize('settings_key', ('fixed_coords', 'FIXED_COORDS'))
def test_fixed_coords(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, calculation, settings_key
):
    """Test the parsing of a successful calculation that has specified a ``fixed_coords`` setting.

    The output files of this test were generated for a calculation of a FCC Si supercell where
    """
    name = f"fixed_coords_{calculation.replace('-', '')}"
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(
        calculation_type=calculation, settings={settings_key: [[True, True, True], [True, True, False]]}
    )
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message


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
    assert orm.Log.collection.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].base.attributes.all,
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
    assert orm.Log.collection.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].base.attributes.all,
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
    assert orm.Log.collection.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_parameters': results['output_parameters'].get_dict(),
        'output_trajectory': results['output_trajectory'].base.attributes.all,
    })


@pytest.mark.parametrize('parameters', (  # yapf:disable
    {'ELECTRONS': {'scf_must_converge': False}},
    {'ELECTRONS': {'electron_maxstep': 0}},
    {'ELECTRONS': {'scf_must_converge': False, 'electron_maxstep': 0}}
))  # yapf:enable
def test_pw_failed_scf_not_converged_intentional(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, parameters
):
    """Test the parsing of an scf calculation that intentionally did not reach convergence.

    This can be the case if `scf_must_converge = False` or `electron_maxstep = 0`.
    """
    name = 'failed_scf_not_converged_intentional'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, generate_inputs(parameters=parameters))
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED.status
    assert orm.Log.collection.get_logs_for(node)
    assert 'output_parameters' in results
    assert 'output_trajectory' in results


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
    assert not orm.Log.collection.get_logs_for(node), [log.message for log in orm.Log.collection.get_logs_for(node)]
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_kpoints': results['output_kpoints'].base.attributes.all,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_structure': results['output_structure'].base.attributes.all,
        'output_trajectory': results['output_trajectory'].base.attributes.all,
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
    assert orm.Log.collection.get_logs_for(node)
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
    assert orm.Log.collection.get_logs_for(node)
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
    assert not orm.Log.collection.get_logs_for(node), [log.message for log in orm.Log.collection.get_logs_for(node)]
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_kpoints': results['output_kpoints'].base.attributes.all,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_structure': results['output_structure'].base.attributes.all,
        'output_trajectory': results['output_trajectory'].base.attributes.all,
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
    assert not orm.Log.collection.get_logs_for(node), [log.message for log in orm.Log.collection.get_logs_for(node)]
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'output_kpoints': results['output_kpoints'].base.attributes.all,
        'output_parameters': results['output_parameters'].get_dict(),
        'output_structure': results['output_structure'].base.attributes.all,
        'output_trajectory': results['output_trajectory'].base.attributes.all,
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
    assert not orm.Log.collection.get_logs_for(node), [log.message for log in orm.Log.collection.get_logs_for(node)]
    assert 'output_parameters' in results
    assert 'output_trajectory' in results
    data_regression.check({
        'energy_vdw': results['output_parameters']['energy_vdw'],
        'array|stress': results['output_trajectory'].base.attributes.all['array|stress'],
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
    assert not orm.Log.collection.get_logs_for(node), [log.message for log in orm.Log.collection.get_logs_for(node)]
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


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_pw_vcrelax_failed_charge_wrong(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename
):
    """Test a `vc-relax` that failed because the integrated charge is different from the expected one."""
    name = f'vcrelax_failed_charge_wrong{filename}'
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
    assert orm.Log.collection.get_logs_for(node)
    assert 'output_parameters' in results


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_pw_vcrelax_failed_symmetry_not_orthogonal(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename
):
    """Test a `vc-relax` that failed because original symmetries no longer map onto new structure."""
    name = f'vcrelax_failed_symmetry_not_orthogonal{filename}'
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
    assert orm.Log.collection.get_logs_for(node)
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
    assert orm.Log.collection.get_logs_for(node)
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
    assert not orm.Log.collection.get_logs_for(node), [log.message for log in orm.Log.collection.get_logs_for(node)]
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
    assert orm.Log.collection.get_logs_for(node)
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
    assert orm.Log.collection.get_logs_for(node)
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
    assert orm.Log.collection.get_logs_for(node)
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
    assert orm.Log.collection.get_logs_for(node)
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
    assert orm.Log.collection.get_logs_for(node)
    assert 'output_kpoints' in results
    assert 'output_parameters' in results
    assert 'output_structure' in results


@pytest.mark.parametrize('filename', ('', '_stdout'))
def test_pw_vcrelax_failed_fft_significant_volume_contraction(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, filename
):
    """Test a `vc-relax` that failed due to significant volume contraction."""
    name = f'vcrelax_failed_fft_significant_volume_contraction{filename}'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    inputs = generate_inputs(calculation_type='vc-relax')
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)
    expected_exit_status = node.process_class.exit_codes.ERROR_RADIAL_FFT_SIGNIFICANT_VOLUME_CONTRACTION.status

    assert calcfunction.is_failed
    assert calcfunction.exit_status == expected_exit_status
    assert orm.Log.collection.get_logs_for(node)
    assert 'output_structure' in results


def test_magnetic_moments_v68(
    fixture_localhost, generate_calc_job_node, generate_parser, generate_inputs, data_regression
):
    """Test the parsing of the magnetic moments in QE v6.8 from stdout."""
    name = 'magnetic_moments_v68'
    entry_point_calc_job = 'quantumespresso.pw'
    entry_point_parser = 'quantumespresso.pw'

    # By setting the `without_xml` option the parsing can only use the stdout content
    inputs = generate_inputs(calculation_type='relax', metadata={'options': {'without_xml': True}})
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, name, inputs)
    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_trajectory' in results

    data_regression.check({
        'atomic_charges':
        results['output_trajectory'].get_array('atomic_charges').tolist(),
        'atomic_magnetic_moments':
        results['output_trajectory'].get_array('atomic_magnetic_moments').tolist(),
    })
