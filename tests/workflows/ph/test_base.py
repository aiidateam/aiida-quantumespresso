# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the `PhBaseWorkChain` class."""
from aiida.common import AttributeDict
from aiida.engine import ProcessHandlerReport

from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.workflows.ph.base import PhBaseWorkChain


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


def test_set_max_seconds(generate_workchain_ph):
    """Test that `max_seconds` gets set in the parameters based on `max_wallclock_seconds` unless already set."""
    inputs = generate_workchain_ph(return_inputs=True)
    max_wallclock_seconds = inputs['ph']['metadata']['options']['max_wallclock_seconds']

    process = generate_workchain_ph(inputs=inputs)
    process.setup()
    process.validate_parameters()
    process.prepare_process()

    expected_max_seconds = max_wallclock_seconds * process.defaults.delta_factor_max_seconds
    assert 'max_seconds' in process.ctx.inputs['parameters']['INPUTPH']
    assert process.ctx.inputs['parameters']['INPUTPH']['max_seconds'] == expected_max_seconds

    # Now check that if `max_seconds` is already explicitly set in the parameters, it is not overwritten.
    inputs = generate_workchain_ph(return_inputs=True)
    max_seconds = 1
    max_wallclock_seconds = inputs['ph']['metadata']['options']['max_wallclock_seconds']
    inputs['ph']['parameters']['INPUTPH']['max_seconds'] = max_seconds

    process = generate_workchain_ph(inputs=inputs)
    process.setup()
    process.validate_parameters()
    process.prepare_process()

    assert 'max_seconds' in process.ctx.inputs['parameters']['INPUTPH']
    assert process.ctx.inputs['parameters']['INPUTPH']['max_seconds'] == max_seconds
