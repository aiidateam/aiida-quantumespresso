"""Tests for the ``NebBaseWorkChain.get_builder_from_protocol`` method."""

import pytest
from aiida.engine import ProcessBuilder

from aiida_quantumespresso.workflows.neb.base import NebBaseWorkChain

pytestmark = pytest.mark.usefixtures('pseudo_family')


def test_get_available_protocols():
    """Test ``NebBaseWorkChain.get_available_protocols``."""
    protocols = NebBaseWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['balanced', 'fast', 'stringent']
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``NebBaseWorkChain.get_default_protocol``."""
    assert NebBaseWorkChain.get_default_protocol() == 'balanced'


def test_default(fixture_code, generate_trajectory):
    """Test ``NebBaseWorkChain.get_builder_from_protocol`` for the default protocol."""
    code = fixture_code('quantumespresso.neb')
    images = generate_trajectory()
    builder = NebBaseWorkChain.get_builder_from_protocol(code, images)

    assert isinstance(builder, ProcessBuilder)
    assert builder.neb.code == code
    assert builder.neb.images == images
    assert sorted(builder.neb.pw.pseudos.keys()) == sorted(images.get_step_structure(-1).get_kind_names())
    assert 'SYSTEM' in builder.neb.pw.parameters.get_dict()
    assert builder.neb.metadata.options['resources']['num_machines'] == 1


def test_overrides_forwarded_to_pw_base(fixture_code, generate_trajectory):
    """Test that ``overrides`` understood by the ``PwBaseWorkChain`` take effect on the builder.

    Note that ``NebBaseWorkChain`` forwards the ``overrides`` to ``PwBaseWorkChain.get_builder_from_protocol`` and then
    transfers a fixed selection of the resulting inputs onto its own builder. Consequently only overrides that are
    understood by the ``PwBaseWorkChain`` spec take effect here; ``NebBaseWorkChain`` specific inputs, such as
    ``neb.parameters`` or ``handler_overrides``, cannot currently be set through the ``overrides``.
    """
    code = fixture_code('quantumespresso.neb')
    images = generate_trajectory()
    overrides = {
        'max_iterations': 3,
        'pw': {'parameters': {'SYSTEM': {'ecutwfc': 42.0}}},
    }
    builder = NebBaseWorkChain.get_builder_from_protocol(code, images, overrides=overrides)

    assert builder.max_iterations == 3
    assert builder.neb.pw.parameters['SYSTEM']['ecutwfc'] == 42.0


def test_options(fixture_code, generate_trajectory):
    """Test specifying ``options`` for the ``get_builder_from_protocol()`` method."""
    code = fixture_code('quantumespresso.neb')
    images = generate_trajectory()

    queue_name = 'super-fast'
    options = {'queue_name': queue_name, 'withmpi': False}
    builder = NebBaseWorkChain.get_builder_from_protocol(code, images, options=options)

    assert builder.neb.metadata.options['queue_name'] == queue_name
    assert builder.neb.metadata.options['withmpi'] is False
