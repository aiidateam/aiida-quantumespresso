# -*- coding: utf-8 -*-
"""Tests for the `PwCalculationTools` class."""
from __future__ import absolute_import
import numpy
import pytest

from aiida import orm
from aiida.common.links import LinkType


def test_pw_get_scf_accuracy(fixture_database, fixture_computer_localhost, generate_calc_job_node):
    """Test the `PwCalculationTools.get_scf_accuracy` method."""
    entry_point_name = 'quantumespresso.pw'

    # Missing `output_trajectory` node
    node = generate_calc_job_node(entry_point_name, fixture_computer_localhost)
    with pytest.raises(ValueError):
        node.tools.get_scf_accuracy()

    # Missing `scf_accuracy` array
    node = generate_calc_job_node(entry_point_name, fixture_computer_localhost)
    trajectory = orm.ArrayData()
    trajectory.add_incoming(node, link_type=LinkType.CREATE, link_label='output_trajectory')
    trajectory.store()

    with pytest.raises(ValueError):
        node.tools.get_scf_accuracy()

    # Missing `scf_accuracy_index` array
    node = generate_calc_job_node(entry_point_name, fixture_computer_localhost)
    trajectory = orm.ArrayData()
    trajectory.set_array('scf_accuracy', numpy.array([1, 1, 1, 2, 2, 2, 2, 2]))
    trajectory.add_incoming(node, link_type=LinkType.CREATE, link_label='output_trajectory')
    trajectory.store()

    with pytest.raises(ValueError):
        node.tools.get_scf_accuracy()

    node = generate_calc_job_node(entry_point_name, fixture_computer_localhost)
    trajectory = orm.ArrayData()
    trajectory.set_array('scf_accuracy', numpy.array([1, 1, 1, 2, 2, 2, 2, 2]))
    trajectory.set_array('scf_accuracy_index', numpy.array([0, 3, 8]))
    trajectory.add_incoming(node, link_type=LinkType.CREATE, link_label='output_trajectory')
    trajectory.store()

    # Invalid indices, there are only two frames
    with pytest.raises(IndexError):
        node.tools.get_scf_accuracy(index=2)
    with pytest.raises(IndexError):
        node.tools.get_scf_accuracy(index=-3)

    node.tools.get_scf_accuracy(index=0) == [1, 1, 1]
    node.tools.get_scf_accuracy(index=1) == [2, 2, 2, 2, 2]
    node.tools.get_scf_accuracy(index=-1) == [2, 2, 2, 2, 2]
    node.tools.get_scf_accuracy(index=-2) == [1, 1, 1]
