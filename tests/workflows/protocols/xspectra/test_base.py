# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the ``XspectraBaseWorkChain.get_builder_from_protocol`` method."""
import io

from aiida.common import LinkType
from aiida.engine import ProcessBuilder
from aiida.plugins import DataFactory, WorkflowFactory
import pytest

XspectraBaseWorkChain = WorkflowFactory('quantumespresso.xspectra.base')


@pytest.fixture
def generate_xspectra_calc_job_node(
    generate_calc_job_node, fixture_localhost, fixture_sandbox, generate_inputs_xspectra
):
    """Generate a ``CalcJobNode`` that would have been created by an ``XspectraCalculation``."""
    from aiida.orm import Dict, FolderData, RemoteData

    def _generate_xspectra_calc_job_node():
        node = generate_calc_job_node(entry_point_name='quantumespresso.xspectra', inputs=generate_inputs_xspectra())

        remote_folder = RemoteData(computer=fixture_localhost, remote_path=fixture_sandbox.abspath)
        remote_folder.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='remote_folder')
        remote_folder.store()

        retrieved = FolderData()
        retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
        retrieved.store()

        output_parameters = Dict()
        output_parameters.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='output_parameters')
        output_parameters.store()

        return node

    return _generate_xspectra_calc_job_node


@pytest.fixture
def get_xspectra_generator_inputs(fixture_code, generate_xspectra_calc_job_node):
    """Generate a set of default inputs for the ``XspectraBaseWorkChain.get_builder_from_protocol()`` method."""
    SinglefileData = DataFactory('core.singlefile')

    parent_calc = generate_xspectra_calc_job_node()

    return {
        'code':
        fixture_code('quantumespresso.xspectra'),
        'parent_folder':
        parent_calc.outputs.remote_folder,
        'core_wfc_data':
        SinglefileData(
            io.StringIO(
                '# number of core states 3 =  1 0;  2 0;'
                '\n6.51344e-05 6.615743462459999e-3'
                '\n6.59537e-05 6.698882211449999e-3'
            )
        ),
        'abs_atom_marker':
        'Si',
    }


def test_get_available_protocols():
    """Test ``XspectraBaseWorkChain.get_available_protocols()``."""
    protocols = XspectraBaseWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['balanced', 'fast', 'replot', 'stringent']
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``XspectraBaseWorkChain.get_default_protocol``."""
    assert XspectraBaseWorkChain.get_default_protocol() == 'balanced'


def test_default(get_xspectra_generator_inputs, data_regression, serialize_builder):
    """Test ``XspectraBaseWorkChain.get_builder_from_protocol`` for the default protocol."""
    # workchain_inputs = get_xspectra_generator_inputs
    builder = XspectraBaseWorkChain.get_builder_from_protocol(**get_xspectra_generator_inputs)

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))


def test_options(get_xspectra_generator_inputs):
    """Test the ``options`` input for the XspectraBaseWorkChain.get_builder_from_protocol() method."""
    queue_name = 'super-fast'
    withmpi = False  # The protocol default is ``True``

    options = {'queue_name': queue_name, 'withmpi': withmpi}
    builder = XspectraBaseWorkChain.get_builder_from_protocol(**get_xspectra_generator_inputs, options=options)

    for subspace in (builder.xspectra.metadata,):
        assert subspace['options']['queue_name'] == queue_name, subspace
