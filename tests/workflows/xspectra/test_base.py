# -*- coding: utf-8 -*-
"""Tests for the XspectraBaseWorkChain."""

from aiida.common import AttributeDict
from aiida.engine import ProcessHandlerReport

from aiida_quantumespresso.calculations.xspectra import XspectraCalculation
from aiida_quantumespresso.workflows.xspectra.base import XspectraBaseWorkChain


def test_setup(generate_workchain_xspectra):
    """Test `XspectraBaseWorkChain.setup`."""
    process = generate_workchain_xspectra()
    process.setup()

    assert isinstance(process.ctx.inputs, AttributeDict)


def test_handle_unrecoverable_failure(generate_workchain_xspectra):
    """Test `XspectraBaseWorkChain.handle_unrecoverable_failure`."""
    process = generate_workchain_xspectra(exit_code=XspectraCalculation.exit_codes.ERROR_NO_RETRIEVED_FOLDER)
    process.setup()

    result = process.handle_unrecoverable_failure(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert result.exit_code == XspectraBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE

    result = process.inspect_process()
    assert result == XspectraBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE


def test_handle_out_of_walltime(generate_workchain_xspectra, fixture_localhost, generate_remote_data):
    """Test `XspectraBaseWorkChain.handle_out_of_walltime`."""
    remote_data = generate_remote_data(computer=fixture_localhost, remote_path='/path/to/remote')
    process = generate_workchain_xspectra(
        exit_code=XspectraCalculation.exit_codes.ERROR_OUT_OF_WALLTIME, xspectra_outputs={'remote_folder': remote_data}
    )
    process.setup()

    result = process.handle_out_of_walltime(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break

    result = process.inspect_process()
    assert result.status == 0


def test_set_time_limit(generate_workchain_xspectra):
    """Test that `time_limit` gets set in the parameters based on `max_wallclock_seconds` unless already set."""
    inputs = generate_workchain_xspectra(return_inputs=True)
    max_wallclock_seconds = inputs['xspectra']['metadata']['options']['max_wallclock_seconds']

    process = generate_workchain_xspectra(inputs=inputs)
    process.setup()
    process.prepare_process()

    expected_time_limit = max_wallclock_seconds * process.defaults.delta_factor_time_limit
    assert 'time_limit' in process.ctx.inputs['parameters']['INPUT_XSPECTRA']
    assert process.ctx.inputs['parameters']['INPUT_XSPECTRA']['time_limit'] == expected_time_limit

    # Now check that if `time_limit` is already explicitly set in the parameters, it is not overwritten.
    inputs = generate_workchain_xspectra(return_inputs=True)
    time_limit = 1
    max_wallclock_seconds = inputs['xspectra']['metadata']['options']['max_wallclock_seconds']
    inputs['xspectra']['parameters']['INPUT_XSPECTRA']['time_limit'] = time_limit

    process = generate_workchain_xspectra(inputs=inputs)
    process.setup()
    process.prepare_process()

    assert 'time_limit' in process.ctx.inputs['parameters']['INPUT_XSPECTRA']
    assert process.ctx.inputs['parameters']['INPUT_XSPECTRA']['time_limit'] == time_limit
