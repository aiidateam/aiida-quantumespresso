# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for the `PdosWorkChain` class."""
from __future__ import absolute_import

from aiida import orm, plugins, engine
from aiida.common import LinkType
from aiida.engine.utils import instantiate_process
from aiida.manage.manager import get_manager

from aiida_quantumespresso.calculations.helpers import pw_input_helper
from aiida_quantumespresso.utils.resources import get_default_options
from aiida_quantumespresso.workflows.pdos.base import PdosWorkChain


def clear_spec():
    # TODO this is required while awaiting fix for aiidateam/aiida-core#3143 # pylint disable=fixme
    if hasattr(PdosWorkChain, '_spec'):
        del PdosWorkChain._spec


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
    fixture_database, generate_code_localhost, fixture_computer_localhost, generate_structure, generate_kpoints_mesh,
    generate_upf_data, generate_calc_job, fixture_sandbox_folder
):
    """Test instantiating the WorkChain, then mock its process, by calling methods in the ``spec.outline``."""
    clear_spec()
    wc_builder = PdosWorkChain.get_builder()

    scf_parameters = {'CONTROL': {'calculation': 'scf'}, 'SYSTEM': {'ecutrho': 240.0, 'ecutwfc': 30.0}}

    scf_inputs = {
        'kpoints': generate_kpoints_mesh(2),
        'pw': {
            'code': generate_code_localhost('quantumespresso.pw', fixture_computer_localhost),
            'structure': generate_structure('Si'),
            'parameters': orm.Dict(dict=scf_parameters),
            'pseudos': {
                'Si': generate_upf_data('Si')
            },
            'metadata': {
                'options': get_default_options()
            }
        }
    }

    dos_parameters = {'Emin': -10, 'Emax': 10, 'DeltaE': 0.01, 'align_to_fermi': True}

    wc_builder.base = scf_inputs
    wc_builder.nscf = {'kpoints': generate_kpoints_mesh(4), 'pw': {'metadata': {'options': get_default_options()}}}
    wc_builder.parameters = dos_parameters
    wc_builder.dos = {
        'code': generate_code_localhost('quantumespresso.dos', fixture_computer_localhost),
        'metadata': {
            'options': get_default_options()
        }
    }
    wc_builder.projwfc = {
        'code': generate_code_localhost('quantumespresso.projwfc', fixture_computer_localhost),
        'metadata': {
            'options': get_default_options()
        }
    }
    wc_builder.test_run = True

    wkchain = instantiate_process_builder(wc_builder)
    assert wkchain.validate() is None

    # run scf
    scf_inputs = wkchain.run_scf()
    scf_wkchain = instantiate_process_cls(plugins.WorkflowFactory('quantumespresso.pw.base'), scf_inputs)
    pw_input_helper(scf_wkchain.inputs.pw.parameters.get_dict(), scf_wkchain.inputs.pw.structure)
    wkchain.ctx.scf_parent_folder = orm.RemoteData(remote_path='/path/on/remote', computer=fixture_computer_localhost)

    # assert wkchain.inspect_scf() is None

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
    remote = orm.RemoteData(remote_path='/path/on/remote', computer=fixture_computer_localhost)
    remote.store()
    remote.add_incoming(mock_wknode, link_type=LinkType.RETURN, link_label='remote_folder')
    result = orm.Dict(dict={'fermi_energy': 6.9029595890428})
    result.store()
    result.add_incoming(mock_wknode, link_type=LinkType.RETURN, link_label='output_parameters')
    wkchain.ctx.workchain_nscf = mock_wknode

    assert wkchain.inspect_nscf() is None

    # mock run dos and projwfc, and check that their inputs are acceptable
    dos_inputs, projwfc_inputs = wkchain.run_dos()
    generate_calc_job(fixture_sandbox_folder, 'quantumespresso.dos', dos_inputs)
    generate_calc_job(fixture_sandbox_folder, 'quantumespresso.projwfc', projwfc_inputs)

    # mock dos & projwfc outputs
    for calc_type in ['dos', 'projwfc']:
        mock_calc = orm.CalcJobNode(computer=fixture_computer_localhost, process_type='quantumespresso.' + calc_type)
        mock_calc.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        mock_calc.set_exit_status(0)
        mock_calc.set_process_state(engine.ProcessState.FINISHED)
        mock_calc.store()
        result = orm.Dict()
        result.add_incoming(mock_calc, link_type=LinkType.CREATE, link_label='output_parameters')
        result.store()
        wkchain.ctx['calc_' + calc_type] = mock_calc

    assert wkchain.inspect_dos() is None

    # store results
    wkchain.results()

    wkchain.update_outputs()
    assert set(wkchain.node.get_outgoing().all_link_labels()) == {
        'projwfc__output_parameters', 'dos__output_parameters', 'nscf__remote_folder', 'nscf__output_parameters'
    }
