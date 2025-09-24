# -*- coding: utf-8 -*-
"""General tests for workflows."""
from aiida.engine import CalcJob
from plumpy import ProcessState
import pytest


@pytest.mark.parametrize(
    'entry_point', (
        'quantumespresso.bands.base',
        'quantumespresso.matdyn.base',
        'quantumespresso.pw.base',
        'quantumespresso.q2r.base',
    )
)
def test_base_unrecoverable_failure(generate_workchain, generate_calc_job_node, generate_inputs, entry_point):
    """Test that the `BaseRestartWorkChain` workflows handle generic unrecoverable failures properly."""
    exit_code = CalcJob.exit_codes.ERROR_SCHEDULER_NODE_FAILURE

    calcjob_entry_point = entry_point.removesuffix('.base')
    calcjob_namespace = calcjob_entry_point.split('.')[-1]
    process = generate_workchain(entry_point, {calcjob_namespace: generate_inputs(calcjob_entry_point)})
    process.setup()

    node = generate_calc_job_node()
    node.set_process_state(ProcessState.FINISHED)
    node.set_exit_status(exit_code.status)

    process.ctx.children = [node]
    process.ctx.iteration = 1

    result = process.inspect_process()
    assert result is None, 'The first inspection should return `None`, i.e. the `BaseRestartWorkChain` should continue'

    node = generate_calc_job_node()
    node.set_process_state(ProcessState.FINISHED)
    node.set_exit_status(exit_code.status)

    process.ctx.children.append(node)
    process.ctx.iteration = 2

    assert process.inspect_process() == process.exit_codes.ERROR_SECOND_CONSECUTIVE_UNHANDLED_FAILURE, \
        'The second inspection should return the proper exit code, i.e. the `BaseRestartWorkChain` should stop'
