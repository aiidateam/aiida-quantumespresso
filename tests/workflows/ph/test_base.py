# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the `PhBaseWorkChain` class."""
from aiida import orm
from aiida.common import AttributeDict, LinkType
from aiida.engine import ProcessHandlerReport, WorkChain
import pytest

from aiida_quantumespresso.calculations.ph import PhCalculation
from aiida_quantumespresso.workflows.ph.base import PhBaseWorkChain


def generate_inputs():
    """Return only those inputs that the parser will expect to be there."""
    return {'parameters': orm.Dict({'INPUTPH': {}})}


@pytest.fixture
def generate_ph_calc_job_node(generate_calc_job_node, fixture_localhost):
    """Generate a ``CalcJobNode`` that would have been created by a ``PhCalculation``."""

    def _generate_ph_calc_job_node():
        node = generate_calc_job_node()

        remote_folder = orm.RemoteData(computer=fixture_localhost, remote_path='/tmp')
        remote_folder.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='remote_folder')
        remote_folder.store()

        retrieved = orm.FolderData()
        retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
        retrieved.store()

        output_parameters = orm.Dict()
        output_parameters.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='output_parameters')
        output_parameters.store()

        return node

    return _generate_ph_calc_job_node


@pytest.mark.usefixtures('aiida_profile')
def test_invalid_inputs(generate_workchain_ph, generate_inputs_ph):
    """Test `PhBaseWorkChain` validation methods."""
    inputs = {'ph': generate_inputs_ph()}
    message = r'Neither `qpoints` nor `qpoints_distance` were specified.'
    with pytest.raises(ValueError, match=message):
        generate_workchain_ph(inputs=inputs)


def test_setup(generate_workchain_ph):
    """Test `PhBaseWorkChain.setup`."""
    process = generate_workchain_ph()
    process.setup()

    assert process.ctx.restart_calc is None
    assert isinstance(process.ctx.inputs, AttributeDict)


@pytest.mark.parametrize(
    ('with_output_structure', 'with_qpoints_distance'),
    ((False, False), (False, True), (True, True)),
)
def test_set_qpoints(generate_workchain_ph, generate_inputs_ph, with_output_structure, with_qpoints_distance):
    """Test `PhBaseWorkChain.set_qpoints`."""
    inputs = {'ph': generate_inputs_ph(with_output_structure=with_output_structure)}
    inputs['qpoints'] = inputs['ph'].pop('qpoints')

    if with_qpoints_distance:
        inputs.pop('qpoints')
        inputs['qpoints_distance'] = orm.Float(0.5)

    process = generate_workchain_ph(inputs=inputs)
    process.setup()
    process.set_qpoints()

    assert 'qpoints' in process.ctx.inputs
    assert isinstance(process.ctx.inputs['qpoints'], orm.KpointsData)

    if not with_qpoints_distance:
        assert process.ctx.inputs['qpoints'] == inputs['qpoints']


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

    assert sorted(process.node.base.links.get_outgoing().all_link_labels()
                  ) == ['output_parameters', 'remote_folder', 'retrieved']


@pytest.mark.parametrize('name', ('merge_outputs', 'merge_outputs_singleq'))
def test_merge_outputs(
    generate_calc_job_node,
    fixture_localhost,
    generate_parser,
    generate_workchain_ph,
    data_regression,
    name,
):
    """Test the ``create_merged_outputs`` step."""

    entry_point_calc_job = 'quantumespresso.ph'
    parser = generate_parser('quantumespresso.ph')
    inputs = generate_inputs()

    node_1 = generate_calc_job_node(
        entry_point_name=entry_point_calc_job, computer=fixture_localhost, test_name=f'{name}_1', inputs=inputs
    )
    results_1, calcjob_1 = parser.parse_from_node(node_1, store_provenance=False)

    results_1['output_parameters'].base.links.add_incoming(
        node_1, link_type=LinkType.CREATE, link_label='output_parameters'
    )
    results_1['output_parameters'].store()

    assert calcjob_1.is_finished, calcjob_1.exception
    assert calcjob_1.exit_status == PhCalculation.exit_codes.ERROR_OUT_OF_WALLTIME.status

    node_2 = generate_calc_job_node(
        entry_point_name=entry_point_calc_job, computer=fixture_localhost, test_name=f'{name}_2', inputs=inputs
    )
    results_2, calcjob_2 = parser.parse_from_node(node_2, store_provenance=False)

    results_2['output_parameters'].base.links.add_incoming(
        node_2, link_type=LinkType.CREATE, link_label='output_parameters'
    )
    results_2['output_parameters'].store()

    assert calcjob_2.is_finished, calcjob_2.exception
    assert calcjob_2.exit_status == 0

    ph_base = generate_workchain_ph()
    ph_base.ctx.children = [node_1, node_2]
    ph_base.ctx.inputs = AttributeDict()
    ph_base.ctx.inputs.parameters = {'INPUTPH': {}}

    ph_base.create_merged_output()
    outputs = ph_base.get_outputs(node_2)

    data_regression.check(outputs['output_parameters'].get_dict())


def test_validate_inputs_excluded_qpoints_distance(generate_inputs_ph):
    """Test that q-points validation passes in case the ports are excluded in a parent work chain."""
    from aiida.engine.utils import instantiate_process
    from aiida.manage.manager import get_manager

    class WrapPhBaseWorkChain(WorkChain):
        """Minimal work chain that wraps a ``PhBaseWorkChain`` for testing excluding q-points inputs."""

        @classmethod
        def define(cls, spec):
            super().define(spec)
            spec.expose_inputs(PhBaseWorkChain, exclude=('qpoints', 'qpoints_distance'))

    inputs = {'ph': generate_inputs_ph()}
    inputs['ph'].pop('qpoints', None)
    runner = get_manager().get_runner()
    instantiate_process(runner, WrapPhBaseWorkChain, **inputs)
