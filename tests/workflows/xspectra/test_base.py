# -*- coding: utf-8 -*-
"""Tests for the XspectraBaseWorkChain.

Currently ToDo:
* Compile required fixtures:
    * Core WFC Data (Done)
* Complete generate_inputs_xspectra (Done)
* Complete generate_workchain_xspectra (Done)
* Finish test_base.py:
    *
"""

from __future__ import absolute_import

from aiida import engine, orm, plugins
from aiida.common import LinkType
from aiida.engine.utils import instantiate_process
from aiida.manage.manager import get_manager
from plumpy import ProcessState

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

def test_default( # pylint: disable=too-many-statements
    generate_workchain_xspectra,
    fixture_localhost,
    generate_remote_data,
    generate_xy_data,
    generate_calc_job,
    generate_calc_job_node,
    fixture_sandbox,
):
    """Test instantiating the main WorkChain and sub-WorkChains, then mock each process in the process spec."""

    xs_workchain = generate_workchain_xspectra()

    assert xs_workchain.setup() is None

    xs_workchain.validate_kpoints()
    assert xs_workchain.validate_kpoints() is None

    # run scf
    scf_inputs = xs_workchain.run_scf()

    # create mock scf_workchain node and outputs
    scf_workchain = instantiate_process_cls(plugins.WorkflowFactory('quantumespresso.pw.base'), scf_inputs)
    pw_input_helper(scf_workchain.inputs.pw.parameters.get_dict(), scf_workchain.inputs.pw.structure)

    scf_workchain_node = scf_workchain.node
    scf_workchain_node.set_exit_status(0)
    scf_workchain_node.set_process_state(ProcessState.FINISHED)
    scf_workchain_node.store()

    remote = generate_remote_data(fixture_localhost, '/path/on/remote', 'quantumespresso.pw')
    remote.store()
    remote.base.links.add_incoming(scf_workchain_node, link_type=LinkType.RETURN, link_label='remote_folder')

    outputs = orm.Dict()
    outputs.store()
    outputs.base.links.add_incoming(scf_workchain_node, link_type=LinkType.RETURN, link_label='output_parameters')

    xs_workchain.ctx.scf_workchain = scf_workchain_node
    # xs_workchain.ctx.scf_parent_folder = remote

    assert xs_workchain.inspect_scf() is None

    # run xspectra production job
    xspectra_prod_inputs = xs_workchain.run_all_xspectra_prod()
    xspectra_prod_params = xspectra_prod_inputs.parameters.get_dict()
    eps_vectors = [
        xspectra_prod_params['INPUT_XSPECTRA']['xepsilon(1)'],
        xspectra_prod_params['INPUT_XSPECTRA']['xepsilon(2)'],
        xspectra_prod_params['INPUT_XSPECTRA']['xepsilon(3)'],
    ]
    # generate_calc_job(fixture_sandbox, 'quantumespresso.xspectra', xspectra_prod_inputs)
    calc_label = 'xas_0'
    prod_label = calc_label + '_prod'
    plot_label = calc_label + '_plot'
    xs_workchain.ctx.xspectra_calc_labels = [calc_label]

    # create mock xspectra production outputs
    entry_point_xs = 'quantumespresso.xspectra'
    mock_calcjob_prod = generate_calc_job_node(
        computer=fixture_localhost, entry_point_name=entry_point_xs, inputs=xspectra_prod_inputs, test_name='default'
    )
    mock_calcjob_prod.label = prod_label + '_iter_1'
    mock_calcjob_prod.set_exit_status(0)
    mock_calcjob_prod.set_process_state(engine.ProcessState.FINISHED)

    prod_out_params = orm.Dict(dict={
        'epsilon_vector': eps_vectors,
        'vector_coord_system': 'crystal',
    })
    prod_out_params.base.links.add_incoming(
        mock_calcjob_prod, link_type=LinkType.CREATE, link_label='output_parameters'
    )
    prod_out_params.store()

    xs_workchain.ctx[prod_label] = mock_calcjob_prod
    # xs_workchain.ctx.finished_lanczos = [mock_calcjob_prod]

    assert xs_workchain.inspect_all_xspectra_prod() is None

    # run xspectra re-plot job
    xspectra_plot_inputs = xs_workchain.run_all_xspectra_plot()
    generate_calc_job(fixture_sandbox, 'quantumespresso.xspectra', xspectra_plot_inputs)

    # create mock xspectra plot outputs
    entry_point_xs = 'quantumespresso.xspectra'
    mock_calcjob_plot = generate_calc_job_node(
        entry_point_name=entry_point_xs, computer=fixture_localhost, inputs=xspectra_plot_inputs
    )
    mock_calcjob_plot.label = plot_label
    mock_calcjob_plot.set_exit_status(0)
    mock_calcjob_plot.set_process_state(engine.ProcessState.FINISHED)

    plot_out_params = orm.Dict(dict={
        'epsilon_vector': eps_vectors,
        'vector_coord_system': 'crystal',
    })
    plot_out_params.base.links.add_incoming(
        mock_calcjob_plot, link_type=LinkType.CREATE, link_label='output_parameters'
    )
    plot_out_params.store()
    plot_out_spectra = generate_xy_data()
    plot_out_spectra.base.links.add_incoming(mock_calcjob_plot, link_type=LinkType.CREATE, link_label='spectra')
    plot_out_spectra.store()

    xs_workchain.ctx[f'{plot_label}'] = mock_calcjob_plot

    assert xs_workchain.inspect_all_xspectra_plot() is None

    xs_workchain.results()

    xs_workchain.update_outputs()

    assert set(xs_workchain.node.base.links.get_outgoing().all_link_labels()
               ) == {'output_parameters_xspectra__xas_0_prod', 'output_parameters_scf', 'output_spectra'}
