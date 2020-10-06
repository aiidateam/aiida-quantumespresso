# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the `PhBaseWorkChain` class."""
import pytest

from plumpy import ProcessState

from aiida.common import AttributeDict
from aiida.engine import ProcessHandlerReport

from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.workflows.ph.base import PhBaseWorkChain


@pytest.fixture
def generate_workchain_ph(generate_workchain, generate_inputs_ph, generate_calc_job_node):
    """Generate an instance of a `PhBaseWorkChain`."""

    def _generate_workchain_ph(exit_code=None):
        entry_point = 'quantumespresso.ph.base'
        process = generate_workchain(entry_point, {'ph': generate_inputs_ph()})

        if exit_code is not None:
            node = generate_calc_job_node()
            node.set_process_state(ProcessState.FINISHED)
            node.set_exit_status(exit_code.status)

            process.ctx.iteration = 1
            process.ctx.children = [node]

        return process

    return _generate_workchain_ph


def test_setup(generate_workchain_ph):
    """Test `PhBaseWorkChain.setup`."""
    process = generate_workchain_ph()
    process.setup()

    assert process.ctx.restart_calc is None
    assert isinstance(process.ctx.inputs, AttributeDict)


def test_handle_unrecoverable_failure(generate_workchain_ph):
    """Test `PhBaseWorkChain.handle_unrecoverable_failure`."""
    process = generate_workchain_ph(exit_code=PhCalculation.exit_codes.ERROR_NO_RETRIEVED_FOLDER)
    process.setup()

    result = process.handle_unrecoverable_failure(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == PhBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE

    result = process.inspect_process()
    assert result == PhBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE


def test_handle_out_of_walltime(generate_workchain_ph):
    """Test `PhBaseWorkChain.handle_out_of_walltime`."""
    process = generate_workchain_ph(exit_code=PhCalculation.exit_codes.ERROR_OUT_OF_WALLTIME)
    process.setup()

    result = process.handle_out_of_walltime(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break

    result = process.inspect_process()
    assert result.status == 0


def test_handle_convergence_not_achieved(generate_workchain_ph):
    """Test `PhBaseWorkChain.handle_convergence_not_achieved`."""
    process = generate_workchain_ph(exit_code=PhCalculation.exit_codes.ERROR_CONVERGENCE_NOT_REACHED)
    process.setup()
    process.validate_parameters()

    alpha_new = PhBaseWorkChain.defaults.alpha_mix * PhBaseWorkChain.defaults.delta_factor_alpha_mix

    result = process.handle_convergence_not_achieved(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert process.ctx.inputs.parameters['INPUTPH']['alpha_mix(1)'] == alpha_new

    result = process.inspect_process()
    assert result.status == 0
