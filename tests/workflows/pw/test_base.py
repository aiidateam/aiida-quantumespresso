# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the `PwBaseWorkChain` class."""
from aiida.common import AttributeDict
from aiida.engine import ExitCode, ProcessHandlerReport
from aiida.orm import Dict
import pytest

from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain


def test_setup(generate_workchain_pw):
    """Test `PwBaseWorkChain.setup`."""
    process = generate_workchain_pw()
    process.setup()

    assert isinstance(process.ctx.inputs, AttributeDict)


def test_handle_unrecoverable_failure(generate_workchain_pw):
    """Test `PwBaseWorkChain.handle_unrecoverable_failure`."""
    process = generate_workchain_pw(exit_code=PwCalculation.exit_codes.ERROR_NO_RETRIEVED_FOLDER)
    process.setup()

    result = process.handle_unrecoverable_failure(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == PwBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE

    result = process.inspect_process()
    assert result == PwBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE


@pytest.mark.parametrize('structure_changed', (
    True,
    False,
))
def test_handle_out_of_walltime(
    generate_workchain_pw, fixture_localhost, generate_remote_data, generate_structure, structure_changed
):
    """Test `PwBaseWorkChain.handle_out_of_walltime`."""
    generate_inputs = {
        'exit_code': PwCalculation.exit_codes.ERROR_OUT_OF_WALLTIME,
        'pw_outputs': {
            'remote_folder': generate_remote_data(computer=fixture_localhost, remote_path='/path/to/remote')
        }
    }
    if structure_changed:
        output_structure = generate_structure()
        generate_inputs['pw_outputs']['output_structure'] = output_structure

    process = generate_workchain_pw(**generate_inputs)
    process.setup()

    result = process.handle_electronic_convergence_not_reached(process.ctx.children[-1])
    result = process.handle_out_of_walltime(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.parameters['CONTROL']['restart_mode'] == 'restart'
    assert result.do_break

    result = process.inspect_process()
    assert result.status == 0

    if structure_changed:
        assert process.ctx.inputs.structure == output_structure


def test_handle_electronic_convergence_not_reached(generate_workchain_pw, fixture_localhost, generate_remote_data):
    """Test `PwBaseWorkChain.handle_electronic_convergence_not_reached`."""
    remote_data = generate_remote_data(computer=fixture_localhost, remote_path='/path/to/remote')

    process = generate_workchain_pw(
        exit_code=PwCalculation.exit_codes.ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED,
        pw_outputs={'remote_folder': remote_data}
    )
    process.setup()

    process.ctx.inputs.parameters['ELECTRONS']['mixing_beta'] = 0.5

    result = process.handle_electronic_convergence_not_reached(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.parameters['ELECTRONS']['mixing_beta'] == \
        process.defaults.delta_factor_mixing_beta * 0.5
    assert process.ctx.inputs.parameters['CONTROL']['restart_mode'] == 'restart'
    assert result.do_break

    result = process.inspect_process()
    assert result.status == 0


@pytest.mark.skip('Reactivate once we have an unrecoverable failure once again.')
def test_handle_known_unrecoverable_failure(generate_workchain_pw):
    """Test `PwBaseWorkChain.handle_known_unrecoverable_failure`."""
    process = generate_workchain_pw()
    process.setup()

    result = process.handle_known_unrecoverable_failure(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == PwBaseWorkChain.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE

    result = process.inspect_process()
    assert result == PwBaseWorkChain.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE


@pytest.mark.parametrize(
    'exit_code', (
        PwCalculation.exit_codes.ERROR_COMPUTING_CHOLESKY,
        PwCalculation.exit_codes.ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED,
        PwCalculation.exit_codes.ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE,
        PwCalculation.exit_codes.ERROR_ZHEGVD_FAILED,
        PwCalculation.exit_codes.ERROR_QR_FAILED,
        PwCalculation.exit_codes.ERROR_EIGENVECTOR_CONVERGENCE,
        PwCalculation.exit_codes.ERROR_BROYDEN_FACTORIZATION,
    )
)
def test_handle_diagonalization_errors(generate_workchain_pw, exit_code):
    """Test `PwBaseWorkChain.handle_diagonalization_errors`."""
    process = generate_workchain_pw(exit_code=exit_code)
    process.setup()

    process.ctx.inputs.parameters['ELECTRONS']['diagonalization'] = 'david'

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.parameters['ELECTRONS']['diagonalization'] == 'ppcg'
    assert result.do_break

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.parameters['ELECTRONS']['diagonalization'] == 'paro'
    assert result.do_break

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.parameters['ELECTRONS']['diagonalization'] == 'cg'
    assert result.do_break

    result = process.inspect_process()
    assert result == PwBaseWorkChain.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE


@pytest.mark.parametrize(
    'exit_code', (
        PwCalculation.exit_codes.ERROR_COMPUTING_CHOLESKY,
        PwCalculation.exit_codes.ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED,
        PwCalculation.exit_codes.ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE,
        PwCalculation.exit_codes.ERROR_ZHEGVD_FAILED,
        PwCalculation.exit_codes.ERROR_QR_FAILED,
        PwCalculation.exit_codes.ERROR_EIGENVECTOR_CONVERGENCE,
        PwCalculation.exit_codes.ERROR_BROYDEN_FACTORIZATION,
    )
)
def test_handle_diagonalization_errors_not_from_david(generate_workchain_pw, exit_code):
    """Test `PwBaseWorkChain.handle_diagonalization_errors` starting from a different diagonalization."""
    process = generate_workchain_pw(exit_code=exit_code)
    process.setup()

    process.ctx.inputs.parameters['ELECTRONS']['diagonalization'] = 'ppcg'

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.parameters['ELECTRONS']['diagonalization'] == 'david'
    assert result.do_break

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.parameters['ELECTRONS']['diagonalization'] == 'paro'
    assert result.do_break

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.parameters['ELECTRONS']['diagonalization'] == 'cg'
    assert result.do_break

    result = process.inspect_process()
    assert result == PwBaseWorkChain.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE


def test_handle_vcrelax_converged_except_final_scf(generate_workchain_pw):
    """Test `PwBaseWorkChain.handle_vcrelax_converged_except_final_scf`."""
    process = generate_workchain_pw(
        exit_code=PwCalculation.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF
    )
    process.setup()

    calculation = process.ctx.children[-1]
    result = process.handle_vcrelax_converged_except_final_scf(calculation)
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == PwBaseWorkChain.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF
    assert process.node.base.links.get_outgoing().all() is not None

    result = process.inspect_process()
    assert result == PwBaseWorkChain.exit_codes.ERROR_IONIC_CONVERGENCE_REACHED_EXCEPT_IN_FINAL_SCF


@pytest.mark.parametrize(
    'exit_code', (
        PwCalculation.exit_codes.ERROR_IONIC_CONVERGENCE_NOT_REACHED,
        PwCalculation.exit_codes.ERROR_IONIC_CYCLE_EXCEEDED_NSTEP,
        PwCalculation.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE,
        PwCalculation.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE,
    )
)
def test_handle_relax_recoverable_ionic_convergence_error(
    generate_workchain_pw, generate_structure, generate_remote_data, fixture_localhost, exit_code
):
    """Test `PwBaseWorkChain.handle_relax_recoverable_ionic_convergence_error`."""
    structure = generate_structure()
    remote_data = generate_remote_data(computer=fixture_localhost, remote_path='/path/to/remote')
    process = generate_workchain_pw(
        pw_outputs={
            'output_structure': structure,
            'remote_folder': remote_data
        }, exit_code=exit_code
    )
    process.setup()

    process.ctx.inputs.parameters['CONTROL']['calculation'] = 'relax'
    process.ctx.inputs.parameters.setdefault('IONS', {})['ion_dynamics'] = 'bfgs'
    result = process.handle_relax_recoverable_ionic_convergence_error(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code.status == 0
    assert process.ctx.inputs.parameters['CONTROL']['restart_mode'] == 'from_scratch'

    result = process.inspect_process()
    assert result.status == 0


@pytest.mark.parametrize(
    'exit_code', (
        PwCalculation.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_FAILURE,
        PwCalculation.exit_codes.ERROR_IONIC_CYCLE_BFGS_HISTORY_AND_FINAL_SCF_FAILURE,
    )
)
def test_handle_relax_recoverable_ionic_convergence_bfgs_history_error(
    generate_workchain_pw, generate_structure, generate_remote_data, fixture_localhost, exit_code
):
    """Test `PwBaseWorkChain.handle_relax_recoverable_ionic_convergence_bfgs_history_error`."""
    structure = generate_structure()
    remote_data = generate_remote_data(computer=fixture_localhost, remote_path='/path/to/remote')
    process = generate_workchain_pw(
        pw_outputs={
            'output_structure': structure,
            'remote_folder': remote_data
        }, exit_code=exit_code
    )
    process.setup()

    # For `relax`, switch to `damp`
    process.ctx.inputs.parameters['CONTROL']['calculation'] = 'relax'
    process.ctx.inputs.parameters.setdefault('IONS', {})['ion_dynamics'] = 'bfgs'
    result = process.handle_relax_recoverable_ionic_convergence_bfgs_history_error(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code.status == 0
    assert process.ctx.inputs.parameters['IONS']['ion_dynamics'] == 'damp'

    # For `vc-relax`, try changing first the `trust_min_radius`
    process.ctx.inputs.parameters['CONTROL']['calculation'] = 'vc-relax'
    process.ctx.inputs.parameters.setdefault('IONS', {})['ion_dynamics'] = 'bfgs'
    process.ctx.inputs.parameters.setdefault('CELL', {})['cell_dynamics'] = 'bfgs'
    result = process.handle_relax_recoverable_ionic_convergence_bfgs_history_error(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code.status == 0
    assert process.ctx.inputs.parameters['CONTROL']['restart_mode'] == 'from_scratch'
    assert process.ctx.inputs.parameters['IONS']['trust_radius_ini'] == 1.0e-3
    assert process.ctx.inputs.parameters['IONS']['trust_radius_min'] == 1.0e-4

    # Then, try `damp` dynamics as a last resort
    result = process.handle_relax_recoverable_ionic_convergence_bfgs_history_error(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code.status == 0
    assert process.ctx.inputs.parameters['IONS']['ion_dynamics'] == 'damp'
    assert process.ctx.inputs.parameters['CELL']['cell_dynamics'] == 'damp-w'


def test_handle_vcrelax_recoverable_fft_significant_volume_contraction_error(generate_workchain_pw, generate_structure):
    """Test `PwBaseWorkChain.handle_vcrelax_recoverable_fft_significant_volume_contraction_error`."""
    exit_code = PwCalculation.exit_codes.ERROR_RADIAL_FFT_SIGNIFICANT_VOLUME_CONTRACTION
    structure = generate_structure()
    process = generate_workchain_pw(pw_outputs={'output_structure': structure}, exit_code=exit_code)
    process.setup()

    result = process.handle_vcrelax_recoverable_fft_significant_volume_contraction_error(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code.status == 0
    assert process.ctx.inputs.parameters['CONTROL']['restart_mode'] == 'from_scratch'
    assert process.ctx.inputs.parameters['CELL']['cell_factor'] == 4

    result = process.inspect_process()
    assert result.status == 0


def test_handle_electronic_convergence_warning(generate_workchain_pw, generate_structure):
    """Test `PwBaseWorkChain.handle_electronic_convergence_warning`."""
    inputs = generate_workchain_pw(return_inputs=True)
    inputs['pw']['parameters']['scf_maxstep'] = 0
    inputs['pw']['parameters']['scf_must_converge'] = False

    structure = generate_structure()

    process = generate_workchain_pw(
        exit_code=PwBaseWorkChain.exit_codes.WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED,
        inputs=inputs,
        pw_outputs={'output_structure': structure}
    )
    process.setup()

    result = process.handle_electronic_convergence_warning(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == PwBaseWorkChain.exit_codes.WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED
    assert process.node.base.links.get_outgoing().all() is not None

    result = process.inspect_process()
    assert result == PwBaseWorkChain.exit_codes.WARNING_ELECTRONIC_CONVERGENCE_NOT_REACHED


def test_sanity_check_no_bands(generate_workchain_pw):
    """Test that `sanity_check_insufficient_bands` does not except if there is no `output_band`, which is optional."""
    process = generate_workchain_pw(exit_code=ExitCode(0))
    process.setup()

    calculation = process.ctx.children[-1]
    assert process.sanity_check_insufficient_bands(calculation) is None


def test_set_max_seconds(generate_workchain_pw):
    """Test that `max_seconds` gets set in the parameters based on `max_wallclock_seconds` unless already set."""
    inputs = generate_workchain_pw(return_inputs=True)
    max_wallclock_seconds = inputs['pw']['metadata']['options']['max_wallclock_seconds']

    process = generate_workchain_pw(inputs=inputs)
    process.setup()
    process.prepare_process()

    expected_max_seconds = max_wallclock_seconds * process.defaults.delta_factor_max_seconds
    assert 'max_seconds' in process.ctx.inputs['parameters']['CONTROL']
    assert process.ctx.inputs['parameters']['CONTROL']['max_seconds'] == expected_max_seconds

    # Now check that if `max_seconds` is already explicitly set in the parameters, it is not overwritten.
    inputs = generate_workchain_pw(return_inputs=True)
    max_seconds = 1
    max_wallclock_seconds = inputs['pw']['metadata']['options']['max_wallclock_seconds']
    inputs['pw']['parameters']['CONTROL']['max_seconds'] = max_seconds

    process = generate_workchain_pw(inputs=inputs)
    process.setup()
    process.prepare_process()

    assert 'max_seconds' in process.ctx.inputs['parameters']['CONTROL']
    assert process.ctx.inputs['parameters']['CONTROL']['max_seconds'] == max_seconds


@pytest.mark.parametrize('restart_mode, expected', (
    ('restart', 'restart'),
    ('from_scratch', 'from_scratch'),
))
def test_restart_mode(generate_workchain_pw, generate_calc_job_node, restart_mode, expected):
    """Test that the ``CONTROL.restart_mode`` specified by the user is always respected."""
    node = generate_calc_job_node('pw', test_name='default')

    inputs = generate_workchain_pw(return_inputs=True)
    inputs['pw']['parent_folder'] = node.outputs.remote_folder
    inputs['pw']['parameters'] = Dict({'CONTROL': {'restart_mode': restart_mode}})

    if restart_mode == 'restart':
        process = generate_workchain_pw(inputs=inputs)
    else:
        with pytest.warns(UserWarning, match='but no input parameters were'):
            process = generate_workchain_pw(inputs=inputs)

    process.setup()

    assert process.ctx.inputs['parameters']['CONTROL']['restart_mode'] == expected
