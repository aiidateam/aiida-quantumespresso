# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the `PwBaseWorkChain` class."""
import pytest

from plumpy import ProcessState

from aiida.common import AttributeDict
from aiida.engine import ProcessHandlerReport

from aiida_quantumespresso.calculations.pw import PwCalculation
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain


@pytest.fixture
def generate_workchain_pw(generate_workchain, generate_inputs_pw, generate_calc_job_node):
    """Generate an instance of a `PwBaseWorkChain`."""

    def _generate_workchain_pw(exit_code=None):
        entry_point = 'quantumespresso.pw.base'
        inputs = generate_inputs_pw()
        kpoints = inputs.pop('kpoints')
        process = generate_workchain(entry_point, {'pw': inputs, 'kpoints': kpoints})

        if exit_code is not None:
            node = generate_calc_job_node()
            node.set_process_state(ProcessState.FINISHED)
            node.set_exit_status(exit_code.status)

            process.ctx.iteration = 1
            process.ctx.children = [node]

        return process

    return _generate_workchain_pw


def test_setup(aiida_profile, generate_workchain_pw):
    """Test `PwBaseWorkChain.setup`."""
    process = generate_workchain_pw()
    process.setup()

    assert process.ctx.restart_calc is None
    assert isinstance(process.ctx.inputs, AttributeDict)


def test_handle_unrecoverable_failure(aiida_profile, generate_workchain_pw):
    """Test `PwBaseWorkChain.handle_unrecoverable_failure`."""
    process = generate_workchain_pw(exit_code=PwCalculation.exit_codes.ERROR_NO_RETRIEVED_FOLDER)
    process.setup()

    result = process.handle_unrecoverable_failure(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == PwBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE

    result = process.inspect_process()
    assert result == PwBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE


def test_handle_out_of_walltime(aiida_profile, generate_workchain_pw):
    """Test `PwBaseWorkChain.handle_out_of_walltime`."""
    process = generate_workchain_pw(exit_code=PwCalculation.exit_codes.ERROR_OUT_OF_WALLTIME)
    process.setup()

    result = process.handle_out_of_walltime(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break

    result = process.inspect_process()
    assert result.status == 0


def test_handle_known_unrecoverable_failure(aiida_profile, generate_workchain_pw):
    """Test `PwBaseWorkChain.handle_known_unrecoverable_failure`."""
    process = generate_workchain_pw(exit_code=PwCalculation.exit_codes.ERROR_COMPUTING_CHOLESKY)
    process.setup()

    result = process.handle_known_unrecoverable_failure(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == PwBaseWorkChain.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE

    result = process.inspect_process()
    assert result == PwBaseWorkChain.exit_codes.ERROR_KNOWN_UNRECOVERABLE_FAILURE


def test_handle_vcrelax_converged_except_final_scf(aiida_profile, generate_workchain_pw):
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
    assert process.node.get_outgoing().all() is not None
    assert process.ctx.restart_calc is calculation

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
    aiida_profile, generate_workchain_pw, generate_structure, exit_code
):
    """Test `PwBaseWorkChain.handle_relax_recoverable_ionic_convergence_error`."""
    from aiida.common.links import LinkType

    process = generate_workchain_pw(exit_code=exit_code)
    process.setup()

    calculation = process.ctx.children[-1]
    structure = generate_structure().store()
    structure.add_incoming(calculation, link_label='output_structure', link_type=LinkType.CREATE)

    result = process.handle_relax_recoverable_ionic_convergence_error(calculation)
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code.status == 0
    assert process.ctx.restart_calc is None

    result = process.inspect_process()
    assert result.status == 0
