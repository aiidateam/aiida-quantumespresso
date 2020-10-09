# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the `Q2rBaseWorkChain` class."""
import pytest

from plumpy import ProcessState

from aiida.common import AttributeDict
from aiida.engine import ProcessHandlerReport

from aiida_quantumespresso.calculations.q2r import Q2rCalculation
from aiida_quantumespresso.workflows.q2r.base import Q2rBaseWorkChain


@pytest.fixture
def generate_workchain_q2r(generate_workchain, generate_inputs_q2r, generate_calc_job_node):
    """Generate an instance of a `Q2rBaseWorkChain`."""

    def _generate_workchain_q2r(exit_code=None):
        entry_point = 'quantumespresso.q2r.base'
        process = generate_workchain(entry_point, {'q2r': generate_inputs_q2r()})

        if exit_code is not None:
            node = generate_calc_job_node()
            node.set_process_state(ProcessState.FINISHED)
            node.set_exit_status(exit_code.status)

            process.ctx.iteration = 1
            process.ctx.children = [node]

        return process

    return _generate_workchain_q2r


def test_setup(generate_workchain_q2r):
    """Test `Q2rBaseWorkChain.setup`."""
    process = generate_workchain_q2r()
    process.setup()

    assert process.ctx.restart_calc is None
    assert isinstance(process.ctx.inputs, AttributeDict)


def test_handle_unrecoverable_failure(generate_workchain_q2r):
    """Test `Q2rBaseWorkChain.handle_unrecoverable_failure`."""
    process = generate_workchain_q2r(exit_code=Q2rCalculation.exit_codes.ERROR_NO_RETRIEVED_FOLDER)
    process.setup()

    result = process.handle_unrecoverable_failure(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == Q2rBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE

    result = process.inspect_process()
    assert result == Q2rBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE
