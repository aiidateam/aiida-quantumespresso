# -*- coding: utf-8 -*-
"""Tests for the `PwCalculationTools` class."""
import numpy as np
import pytest

from aiida import orm
from aiida.common.links import LinkType


def test_pw_get_scf_accuracy(fixture_localhost, generate_calc_job_node):
    """Test the `PwCalculationTools.get_scf_accuracy` method."""
    entry_point_name = 'quantumespresso.pw'

    # Missing `output_trajectory` node
    node = generate_calc_job_node(entry_point_name, fixture_localhost)
    with pytest.raises(ValueError):
        node.tools.get_scf_accuracy()

    # Missing `scf_accuracy` array
    node = generate_calc_job_node(entry_point_name, fixture_localhost)
    trajectory = orm.ArrayData()
    trajectory.add_incoming(node, link_type=LinkType.CREATE, link_label='output_trajectory')
    trajectory.store()

    with pytest.raises(ValueError):
        node.tools.get_scf_accuracy()

    # Missing `scf_accuracy_index` array
    node = generate_calc_job_node(entry_point_name, fixture_localhost)
    trajectory = orm.ArrayData()
    trajectory.set_array('scf_accuracy', np.array([1, 1, 1, 2, 2, 2, 2, 2]))
    trajectory.add_incoming(node, link_type=LinkType.CREATE, link_label='output_trajectory')
    trajectory.store()

    with pytest.raises(ValueError):
        node.tools.get_scf_accuracy()

    node = generate_calc_job_node(entry_point_name, fixture_localhost)
    trajectory = orm.ArrayData()
    trajectory.set_array('scf_accuracy', np.array([1, 1, 1, 2, 2, 2, 2, 2]))
    trajectory.set_array('scf_iterations', np.array([3, 5]))
    trajectory.add_incoming(node, link_type=LinkType.CREATE, link_label='output_trajectory')
    trajectory.store()

    # Invalid indices, there are only two frames
    with pytest.raises(IndexError):
        node.tools.get_scf_accuracy(index=2)
    with pytest.raises(IndexError):
        node.tools.get_scf_accuracy(index=-3)

    assert np.array_equal(node.tools.get_scf_accuracy(index=0), np.array([1, 1, 1]))
    assert np.array_equal(node.tools.get_scf_accuracy(index=1), np.array([2, 2, 2, 2, 2]))
    assert np.array_equal(node.tools.get_scf_accuracy(index=-1), np.array([2, 2, 2, 2, 2]))
    assert np.array_equal(node.tools.get_scf_accuracy(index=-2), np.array([1, 1, 1]))
