# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the `BandsBaseWorkChain` class."""
from aiida.common import AttributeDict
from aiida.engine import ProcessHandlerReport
from plumpy import ProcessState
import pytest

from aiida_quantumespresso.calculations.bands import BandsCalculation
from aiida_quantumespresso.workflows.bands.base import BandsBaseWorkChain


@pytest.fixture
def generate_workchain_bands(generate_workchain, generate_inputs_bands, generate_calc_job_node):
    """Generate an instance of a `BandsBaseWorkChain`."""

    def _generate_workchain_bands(exit_code=None):
        entry_point = 'quantumespresso.bands.base'
        process = generate_workchain(entry_point, {'bands': generate_inputs_bands()})

        if exit_code is not None:
            node = generate_calc_job_node()
            node.set_process_state(ProcessState.FINISHED)
            node.set_exit_status(exit_code.status)

            process.ctx.iteration = 1
            process.ctx.children = [node]

        return process

    return _generate_workchain_bands


def test_setup(generate_workchain_bands):
    """Test `BandsBaseWorkChain.setup`."""
    process = generate_workchain_bands()
    process.setup()

    assert process.ctx.restart_calc is None
    assert isinstance(process.ctx.inputs, AttributeDict)


def test_handle_unrecoverable_failure(generate_workchain_bands):
    """Test `BandsBaseWorkChain.handle_unrecoverable_failure`."""
    process = generate_workchain_bands(exit_code=BandsCalculation.exit_codes.ERROR_NO_RETRIEVED_FOLDER)
    process.setup()

    result = process.handle_unrecoverable_failure(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == BandsBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE

    result = process.inspect_process()
    assert result == BandsBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE
