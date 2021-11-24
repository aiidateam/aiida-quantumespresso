# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the `PhBaseWorkChain` class."""
import pytest

from aiida.common import AttributeDict, LinkType
from aiida.engine import ProcessHandlerReport
from aiida.orm import RemoteData, FolderData, Dict

from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.workflows.ph.base import PhBaseWorkChain


@pytest.fixture
def generate_ph_calc_job_node(generate_calc_job_node, fixture_localhost):
    """Generate a ``CalcJobNode`` that would have been created by a ``PhCalculation``."""

    def _generate_ph_calc_job_node():
        node = generate_calc_job_node()

        remote_folder = RemoteData(computer=fixture_localhost, remote_path='/tmp')
        remote_folder.add_incoming(node, link_type=LinkType.CREATE, link_label='remote_folder')
        remote_folder.store()

        retrieved = FolderData()
        retrieved.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
        retrieved.store()

        output_parameters = Dict()
        output_parameters.add_incoming(node, link_type=LinkType.CREATE, link_label='output_parameters')
        output_parameters.store()

        return node

    return _generate_ph_calc_job_node


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


def test_handle_scheduler_out_of_walltime(generate_workchain_ph):
    """Test `PhBaseWorkChain.handle_scheduler_out_of_walltime`."""
    inputs = generate_workchain_ph(return_inputs=True)
    max_wallclock_seconds = inputs['ph']['metadata']['options']['max_wallclock_seconds']
    max_seconds = max_wallclock_seconds * PhBaseWorkChain.defaults.delta_factor_max_seconds

    process = generate_workchain_ph(exit_code=PhCalculation.exit_codes.ERROR_SCHEDULER_OUT_OF_WALLTIME)
    process.setup()
    process.validate_parameters()
    process.prepare_process()

    max_seconds_new = max_seconds * 0.5

    result = process.handle_scheduler_out_of_walltime(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert process.ctx.inputs.parameters['INPUTPH']['max_seconds'] == max_seconds_new
    assert not process.ctx.inputs.parameters['INPUTPH']['recover']

    result = process.inspect_process()
    assert result.status == 0


def test_handle_convergence_not_reached(generate_workchain_ph):
    """Test `PhBaseWorkChain.handle_convergence_not_reached`."""
    process = generate_workchain_ph(exit_code=PhCalculation.exit_codes.ERROR_CONVERGENCE_NOT_REACHED)
    process.setup()
    process.validate_parameters()

    alpha_new = PhBaseWorkChain.defaults.alpha_mix * PhBaseWorkChain.defaults.delta_factor_alpha_mix

    result = process.handle_convergence_not_reached(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert result.do_break
    assert process.ctx.inputs.parameters['INPUTPH']['alpha_mix(1)'] == alpha_new

    result = process.inspect_process()
    assert result.status == 0


def test_handle_diagonalization_errors(generate_workchain_ph):
    """Test `PhBaseWorkChain.handle_diagonalization_errors`."""
    process = generate_workchain_ph(exit_code=PhCalculation.exit_codes.ERROR_COMPUTING_CHOLESKY)
    process.setup()
    process.validate_parameters()
    process.prepare_process()

    process.ctx.inputs.parameters['INPUTPH']['diagonalization'] = 'david'

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert isinstance(result, ProcessHandlerReport)
    assert process.ctx.inputs.parameters['INPUTPH']['diagonalization'] == 'cg'
    assert result.do_break

    result = process.handle_diagonalization_errors(process.ctx.children[-1])
    assert result.do_break

    result = process.inspect_process()
    assert result == PhBaseWorkChain.exit_codes.ERROR_UNRECOVERABLE_FAILURE


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


def test_results(generate_workchain_ph, generate_ph_calc_job_node):
    """Test `PhBaseWorkChain.results`."""
    process = generate_workchain_ph(exit_code=PhCalculation.exit_codes.ERROR_CONVERGENCE_NOT_REACHED)
    process.setup()
    process.ctx.children = [generate_ph_calc_job_node()]

    # Now call the ``results`` method that should check the outputs. Need to call ``update_outputs`` to actually have
    # the output links generated.
    process.results()
    process.update_outputs()

    assert sorted(process.node.get_outgoing().all_link_labels()) == ['output_parameters', 'remote_folder', 'retrieved']
