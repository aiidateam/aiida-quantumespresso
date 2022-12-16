# -*- coding: utf-8 -*-
"""Tests for the `PdosWorkChain` class."""
from __future__ import absolute_import

import io

from aiida import engine, orm
from aiida.common import LinkType
from aiida.engine.utils import instantiate_process
from aiida.manage.manager import get_manager
from plumpy import ProcessState
import pytest

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


@pytest.fixture
def generate_workchain_xspectra_core(generate_inputs_pw, generate_workchain, generate_inputs_xspectra):
    """Generate an instance of a `XspectraCoreWorkChain`."""

    def _generate_workchain_xspectra_core():
        from aiida.orm import Bool, List, SinglefileData, Str

        entry_point = 'quantumespresso.xspectra.core'

        scf_pw_inputs = generate_inputs_pw()
        xs_prod_inputs = generate_inputs_xspectra()

        xs_prod = {'xspectra': xs_prod_inputs}
        xs_prod_inputs.pop('parent_folder')
        xs_prod_inputs.pop('kpoints')

        kpoints = scf_pw_inputs.pop('kpoints')
        structure = scf_pw_inputs.pop('structure')
        scf = {'pw': scf_pw_inputs, 'kpoints': kpoints}

        inputs = {
            'structure':
            structure,
            'scf':
            scf,
            'xs_prod':
            xs_prod,
            'xs_plot':
            xs_prod,
            'run_replot':
            Bool(False),
            'dry_run':
            Bool(True),
            'abs_atom_marker':
            Str('Si'),
            'core_wfc_data':
            SinglefileData(
                io.StringIO(
                    '# number of core states 3 =  1 0;  2 0;'
                    '\n6.51344e-05 6.615743462459999e-3'
                    '\n6.59537e-05 6.698882211449999e-3'
                )
            ),
            'eps_vectors':
            List(list=[[1., 0., 0.]])
        }

        return generate_workchain(entry_point, inputs)

    return _generate_workchain_xspectra_core


def test_default(
    generate_inputs_xspectra,
    generate_workchain_pw,
    generate_workchain_xspectra,
    generate_workchain_xspectra_core, #pylint: disable=redefined-outer-name
    fixture_localhost,
    generate_remote_data,
    generate_calc_job_node,
    generate_xy_data,
):
    """Test instantiating the WorkChain, then mock its process, by calling methods in the ``spec.outline``."""

    wkchain = generate_workchain_xspectra_core()

    assert wkchain.setup() is None
    assert wkchain.should_run_upf2plotcore() is False
    assert wkchain.should_run_replot() is False

    # run scf
    scf_inputs = wkchain.run_scf()

    scf_wkchain = generate_workchain_pw(inputs=scf_inputs)
    scf_wkchain.node.set_process_state(ProcessState.FINISHED)
    scf_wkchain.node.set_exit_status(0)
    pw_input_helper(scf_wkchain.inputs.pw.parameters.get_dict(), scf_wkchain.inputs.pw.structure)

    remote = generate_remote_data(computer=fixture_localhost, remote_path='/path/on/remote')
    remote.store()
    remote.base.links.add_incoming(scf_wkchain.node, link_type=LinkType.RETURN, link_label='remote_folder')

    result = orm.Dict()
    result.store()
    result.base.links.add_incoming(scf_wkchain.node, link_type=LinkType.RETURN, link_label='output_parameters')

    wkchain.ctx.scf_workchain = scf_wkchain.node
    wkchain.ctx.scf_parent_folder = remote

    assert wkchain.inspect_scf() is None

    # mock run the xs_prod step
    xs_prod_inputs = wkchain.run_all_xspectra_prod()
    # mock xs_prod outputs
    xs_prod_wc = generate_workchain_xspectra(inputs=xs_prod_inputs)
    xs_prod_node = xs_prod_wc.node
    xs_prod_node.label = 'xas_0_prod'
    xs_prod_node.store()
    xs_prod_node.set_exit_status(0)
    xs_prod_node.set_process_state(engine.ProcessState.FINISHED)

    # mock an XspectraCalculation node, as the WorkChain needs one, since it's looking for a
    # "creator" node in the final steps
    xspectra_node = generate_calc_job_node(
        entry_point_name='quantumespresso.xspectra', inputs=generate_inputs_xspectra()
    )
    xspectra_node.store()

    result = orm.Dict(dict={'xepsilon': [1., 0., 0.], 'xcoordcrys': False})
    result.store()
    result.base.links.add_incoming(xs_prod_node, link_type=LinkType.RETURN, link_label='output_parameters')
    result.base.links.add_incoming(xspectra_node, link_type=LinkType.CREATE, link_label='output_parameters')

    spectra = generate_xy_data()
    spectra.store()
    spectra.base.links.add_incoming(xs_prod_node, link_type=LinkType.RETURN, link_label='spectra')
    spectra.base.links.add_incoming(xspectra_node, link_type=LinkType.CREATE, link_label='spectra')

    wkchain.ctx['xas_0_prod'] = xs_prod_node
    wkchain.ctx.xspectra_calc_labels = ['xas_0']

    assert wkchain.inspect_all_xspectra_prod() is None
    assert wkchain.ctx.all_lanczos_computed is True

    # process results
    wkchain.results()

    wkchain.update_outputs()

    assert set(wkchain.node.base.links.get_outgoing().all_link_labels()
               ) == {'output_parameters_scf', 'output_parameters_xspectra__xas_0_prod', 'output_spectra'}
