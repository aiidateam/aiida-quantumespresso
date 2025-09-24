# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the `PwBaseWorkChain` class."""
from aiida.common import AttributeDict
from aiida.engine import ProcessHandlerReport
import pytest

from aiida_quantumespresso.calculations.neb import NebCalculation
from aiida_quantumespresso.workflows.neb.base import NebBaseWorkChain


def test_setup(generate_workchain_neb):
    """Test `NebBaseWorkChain.setup`."""
    process = generate_workchain_neb()
    process.setup()

    assert isinstance(process.ctx.inputs, AttributeDict)


def test_handle_electronic_convergence_not_reached(generate_workchain_neb, fixture_localhost, generate_remote_data):
    """Test `NebBaseWorkChain.handle_electronic_convergence_not_reached`."""
    remote_data = generate_remote_data(computer=fixture_localhost, remote_path='/path/to/remote')

    process = generate_workchain_neb(
        exit_code=NebCalculation.exit_codes.ERROR_ELECTRONIC_CONVERGENCE_NOT_REACHED,
        neb_outputs={'remote_folder': remote_data}
    )
    process.setup()

    process.ctx.inputs.pw.parameters['ELECTRONS']['mixing_beta'] = 0.5

    result = process.handle_electronic_convergence_not_reached(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.pw.parameters['ELECTRONS']['mixing_beta'] == \
        process.defaults.delta_factor_mixing_beta * 0.5
    assert process.ctx.inputs.parameters['PATH']['restart_mode'] == 'restart'
    assert result.do_break

    result = process.inspect_process()
    assert result.status == 0


@pytest.mark.parametrize(
    'exit_code', (
        NebCalculation.exit_codes.ERROR_COMPUTING_CHOLESKY,
        NebCalculation.exit_codes.ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED,
        NebCalculation.exit_codes.ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE,
        NebCalculation.exit_codes.ERROR_ZHEGVD_FAILED,
        NebCalculation.exit_codes.ERROR_QR_FAILED,
        NebCalculation.exit_codes.ERROR_EIGENVECTOR_CONVERGENCE,
        NebCalculation.exit_codes.ERROR_BROYDEN_FACTORIZATION,
    )
)
def test_handle_diagonalization_errors(generate_workchain_neb, exit_code):
    """Test `NebBaseWorkChain.handle_diagonalization_errors`."""
    process = generate_workchain_neb(exit_code=exit_code)
    process.setup()

    process.ctx.inputs.pw.parameters['ELECTRONS']['diagonalization'] = 'david'

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.pw.parameters['ELECTRONS']['diagonalization'] == 'ppcg'
    assert result.do_break

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.pw.parameters['ELECTRONS']['diagonalization'] == 'paro'
    assert result.do_break

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.pw.parameters['ELECTRONS']['diagonalization'] == 'cg'
    assert result.do_break

    result = process.inspect_process()
    assert result == NebBaseWorkChain.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE


@pytest.mark.parametrize(
    'exit_code', (
        NebCalculation.exit_codes.ERROR_COMPUTING_CHOLESKY,
        NebCalculation.exit_codes.ERROR_DIAGONALIZATION_TOO_MANY_BANDS_NOT_CONVERGED,
        NebCalculation.exit_codes.ERROR_S_MATRIX_NOT_POSITIVE_DEFINITE,
        NebCalculation.exit_codes.ERROR_ZHEGVD_FAILED,
        NebCalculation.exit_codes.ERROR_QR_FAILED,
        NebCalculation.exit_codes.ERROR_EIGENVECTOR_CONVERGENCE,
        NebCalculation.exit_codes.ERROR_BROYDEN_FACTORIZATION,
    )
)
def test_handle_diagonalization_errors_not_from_david(generate_workchain_neb, exit_code):
    """Test `NebBaseWorkChain.handle_diagonalization_errors` starting from a different diagonalization."""
    process = generate_workchain_neb(exit_code=exit_code)
    process.setup()

    process.ctx.inputs.pw.parameters['ELECTRONS']['diagonalization'] = 'ppcg'

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.pw.parameters['ELECTRONS']['diagonalization'] == 'david'
    assert result.do_break

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.pw.parameters['ELECTRONS']['diagonalization'] == 'paro'
    assert result.do_break

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.pw.parameters['ELECTRONS']['diagonalization'] == 'cg'
    assert result.do_break

    result = process.inspect_process()
    assert result == NebBaseWorkChain.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE


@pytest.mark.parametrize('exit_code', (NebCalculation.exit_codes.ERROR_NEB_INTERRUPTED_PARTIAL_TRAJECTORY,))
def test_handle_neb_interrupted_partial_trajectory(
    generate_workchain_neb, generate_remote_data, fixture_localhost, exit_code
):
    """Test `NebBaseWorkChain.handle_neb_interrupted_partial_trajectory`."""
    remote_data = generate_remote_data(computer=fixture_localhost, remote_path='/path/to/remote')
    process = generate_workchain_neb(neb_outputs={'remote_folder': remote_data}, exit_code=exit_code)
    process.setup()

    result = process.handle_neb_interrupted_partial_trajectory(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code.status == 0
    assert process.ctx.inputs.parameters['PATH']['restart_mode'] == 'restart'

    result = process.inspect_process()
    assert result.status == 0


@pytest.mark.parametrize('exit_code', (NebCalculation.exit_codes.ERROR_NEB_CYCLE_EXCEEDED_NSTEP,))
def test_handle_neb_cycle_exceeded_nstep_error(
    generate_workchain_neb, generate_remote_data, fixture_localhost, exit_code
):
    """Test `NebBaseWorkChain.handle_neb_cycle_exceeded_nstep`."""
    remote_data = generate_remote_data(computer=fixture_localhost, remote_path='/path/to/remote')
    process = generate_workchain_neb(neb_outputs={'remote_folder': remote_data}, exit_code=exit_code)
    process.setup()
    input_nsteps = process.ctx.inputs.parameters['PATH']['nstep_path'] if 'nstep_path' in process.ctx.inputs.parameters[
        'PATH'] else 1
    result = process.handle_neb_cycle_exceeded_nstep(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code.status == 0
    assert process.ctx.inputs.parameters['PATH']['restart_mode'] == 'restart'
    assert process.ctx.inputs.parameters['PATH']['nstep_path'] == input_nsteps * 2
    result = process.inspect_process()
    assert result.status == 0
