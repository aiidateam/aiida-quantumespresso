# -*- coding: utf-8 -*-
"""Unit tests for the :py:mod:`~aiida_quantumespresso.utils.restart` module."""
import pytest

from aiida.engine import ProcessBuilder
from aiida_quantumespresso.utils import restart


def generate_inputs():
    """Return a dictionary of inputs for a `CalcJobNode` fixture to be created."""
    from aiida import orm
    return {'parameters': orm.Dict(dict={}), 'settings': orm.Dict(dict={})}


def test_restart(fixture_localhost, generate_calc_job_node):
    """Test the generics of the `get_builder_restart`."""
    entry_point_calc_job = 'quantumespresso.dos'
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default', generate_inputs())

    # The `DosCalculation` class is not yet supported
    with pytest.raises(TypeError):
        restart.get_builder_restart(node)


def test_restart_cp(fixture_localhost, generate_calc_job_node):
    """Test the `get_builder_restart` for a completed `CpCalculation`."""
    entry_point_calc_job = 'quantumespresso.cp'
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default', generate_inputs())

    builder = restart.get_builder_restart(node)
    parameters = builder.parameters.get_dict()

    assert isinstance(builder, ProcessBuilder)
    assert parameters['CONTROL']['restart_mode'] == 'restart'

    # Force `from_scratch`
    builder = restart.get_builder_restart(node, from_scratch=True)
    parameters = builder.parameters.get_dict()

    assert isinstance(builder, ProcessBuilder)
    assert parameters['CONTROL']['restart_mode'] == 'from_scratch'


def test_restart_neb(fixture_localhost, generate_calc_job_node):
    """Test the `get_builder_restart` for a completed `NebCalculation`."""
    entry_point_calc_job = 'quantumespresso.neb'
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default', generate_inputs())

    builder = restart.get_builder_restart(node)
    parameters = builder.parameters.get_dict()

    assert isinstance(builder, ProcessBuilder)
    assert parameters['PATH']['restart_mode'] == 'restart'

    # Force `from_scratch`
    builder = restart.get_builder_restart(node, from_scratch=True)
    parameters = builder.parameters.get_dict()

    assert isinstance(builder, ProcessBuilder)
    assert parameters['PATH']['restart_mode'] == 'from_scratch'


def test_restart_ph(fixture_localhost, generate_calc_job_node):
    """Test the `get_builder_restart` for a completed `PhCalculation`."""
    entry_point_calc_job = 'quantumespresso.ph'
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default', generate_inputs())

    builder = restart.get_builder_restart(node)
    parameters = builder.parameters.get_dict()

    assert isinstance(builder, ProcessBuilder)
    assert parameters['INPUTPH']['recover'] is True


def test_restart_pw(fixture_localhost, generate_calc_job_node):
    """Test the `get_builder_restart` for a completed `PwCalculation`."""
    entry_point_calc_job = 'quantumespresso.pw'
    node = generate_calc_job_node(entry_point_calc_job, fixture_localhost, 'default', generate_inputs())

    builder = restart.get_builder_restart(node)
    parameters = builder.parameters.get_dict()

    assert isinstance(builder, ProcessBuilder)
    assert parameters['CONTROL']['restart_mode'] == 'restart'

    # Force `from_scratch`
    builder = restart.get_builder_restart(node, from_scratch=True)
    parameters = builder.parameters.get_dict()

    assert isinstance(builder, ProcessBuilder)
    assert parameters['CONTROL']['restart_mode'] == 'from_scratch'
