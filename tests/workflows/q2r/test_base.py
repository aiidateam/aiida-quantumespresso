"""Tests for the `Q2rBaseWorkChain` class."""

import pytest
from aiida.common import AttributeDict
from plumpy import ProcessState


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
