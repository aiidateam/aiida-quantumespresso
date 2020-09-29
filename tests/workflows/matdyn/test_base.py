# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the `MatdynBaseWorkChain` class."""
import pytest

from plumpy import ProcessState

from aiida.common import AttributeDict
from aiida.engine import ProcessHandlerReport

from aiida_quantumespresso.calculations.matdyn import MatdynCalculation
from aiida_quantumespresso.workflows.matdyn.base import MatdynBaseWorkChain


@pytest.fixture
def generate_workchain_matdyn(generate_workchain, generate_inputs_matdyn, generate_calc_job_node):
    """Generate an instance of a `MatdynBaseWorkChain`."""

    def _generate_workchain_matdyn(exit_code=None):
        entry_point = 'quantumespresso.matdyn.base'
        process = generate_workchain(entry_point, {'matdyn': generate_inputs_matdyn()})

        if exit_code is not None:
            node = generate_calc_job_node()
            node.set_process_state(ProcessState.FINISHED)
            node.set_exit_status(exit_code.status)

            process.ctx.iteration = 1
            process.ctx.children = [node]

        return process

    return _generate_workchain_matdyn


def test_setup(generate_workchain_matdyn):
    """Test `MatdynBaseWorkChain.setup`."""
    process = generate_workchain_matdyn()
    process.setup()

    assert process.ctx.restart_calc is None
    assert isinstance(process.ctx.inputs, AttributeDict)


def test_handle_unrecoverable_failure(generate_workchain_matdyn):
    """Test `MatdynBaseWorkChain.handle_unrecoverable_failure`."""
    process = generate_workchain_matdyn(exit_code=MatdynCalculation.exit_codes.ERROR_NO_RETRIEVED_FOLDER)
    process.setup()

    result = process.handle_unrecoverable_failure(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == MatdynBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE

    result = process.inspect_process()
    assert result == MatdynBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE
