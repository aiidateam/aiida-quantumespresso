# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `PdosWorkChain` class."""
from __future__ import absolute_import

from plumpy import ProcessState

from aiida import orm, plugins, engine
from aiida.common import LinkType
from aiida.engine.utils import instantiate_process
from aiida.manage.manager import get_manager

from aiida_quantumespresso.calculations.helpers import pw_input_helper


def instantiate_process_builder(builder):
    """Instantiate a process, from a `ProcessBuilder`."""
    manager = get_manager()
    runner = manager.get_runner()
    return instantiate_process(runner, builder)


def instantiate_process_cls(process_cls, inputs):
    """Instantiate a process, from its class and inputs."""
    manager = get_manager()
    runner = manager.get_runner()
    return instantiate_process(runner, process_cls, **inputs)


def test_default(
    generate_workchain_pdos,
    generate_workchain_pw,
    fixture_localhost,
    generate_remote_data,
    generate_calc_job,
    generate_calc_job_node,
    fixture_sandbox,
    generate_bands_data,
):
    """Test instantiating the WorkChain, then mock its process, by calling methods in the ``spec.outline``."""

    wkchain = generate_workchain_pdos()
    assert wkchain.setup() is None
    assert wkchain.serial_clean() is False

    # run scf
    scf_inputs = wkchain.run_scf()

    scf_wkchain = generate_workchain_pw(inputs=scf_inputs)
    scf_wkchain.node.set_process_state(ProcessState.FINISHED)
    scf_wkchain.node.set_exit_status(0)
    pw_input_helper(scf_wkchain.inputs.pw.parameters.get_dict(), scf_wkchain.inputs.pw.structure)

    remote = generate_remote_data(computer=fixture_localhost, remote_path='/path/on/remote')
    remote.store()
    remote.add_incoming(scf_wkchain.node, link_type=LinkType.RETURN, link_label='remote_folder')

    wkchain.ctx.workchain_scf = scf_wkchain.node
    wkchain.ctx.scf_parent_folder = remote

    assert wkchain.inspect_scf() is None

    # run nscf
    nscf_inputs = wkchain.run_nscf()

    # mock nscf outputs
    # TODO ensure this test fails if the output link from PwCalculation changes from `output_parameters` # pylint: disable=fixme
    mock_workchain = instantiate_process_cls(plugins.WorkflowFactory('quantumespresso.pw.base'), nscf_inputs)
    pw_input_helper(mock_workchain.inputs.pw.parameters.get_dict(), mock_workchain.inputs.pw.structure)

    mock_wknode = mock_workchain.node
    mock_wknode.set_exit_status(0)
    mock_wknode.set_process_state(engine.ProcessState.FINISHED)
    mock_wknode.store()

    remote = generate_remote_data(computer=fixture_localhost, remote_path='/path/on/remote')
    remote.store()
    remote.add_incoming(mock_wknode, link_type=LinkType.RETURN, link_label='remote_folder')

    result = orm.Dict(dict={'fermi_energy': 6.9029595890428})
    result.store()
    result.add_incoming(mock_wknode, link_type=LinkType.RETURN, link_label='output_parameters')

    bands_data = generate_bands_data()
    bands_data.store()
    bands_data.add_incoming(mock_wknode, link_type=LinkType.RETURN, link_label='output_band')

    wkchain.ctx.workchain_nscf = mock_wknode

    assert wkchain.inspect_nscf() is None

    # mock run dos and projwfc, and check that their inputs are acceptable
    dos_inputs, projwfc_inputs = wkchain.run_pdos_parallel()
    generate_calc_job(fixture_sandbox, 'quantumespresso.dos', dos_inputs)
    generate_calc_job(fixture_sandbox, 'quantumespresso.projwfc', projwfc_inputs)

    # mock dos & projwfc outputs
    for calc_type in ['dos', 'projwfc']:
        entry_point = 'quantumespresso.' + calc_type
        mock_calc = generate_calc_job_node(entry_point_name=entry_point, computer=fixture_localhost)
        mock_calc.set_exit_status(0)
        mock_calc.set_process_state(engine.ProcessState.FINISHED)

        result = orm.Dict()
        result.add_incoming(mock_calc, link_type=LinkType.CREATE, link_label='output_parameters')
        result.store()

        wkchain.ctx['calc_' + calc_type] = mock_calc

    assert wkchain.inspect_dos_serial() is None
    assert wkchain.inspect_projwfc_serial() is None
    assert wkchain.inspect_pdos_parallel() is None

    # store results
    wkchain.results()

    wkchain.update_outputs()

    assert set(wkchain.node.get_outgoing().all_link_labels()) == {
        'projwfc__output_parameters', 'dos__output_parameters', 'nscf__remote_folder', 'nscf__output_parameters',
        'nscf__output_band'
    }
