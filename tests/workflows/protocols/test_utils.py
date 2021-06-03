# -*- coding: utf-8 -*-
"""Tests for the utility functions for the protocols."""


def test_recursive_merge():
    """Test the recursive merge function."""
    from aiida_quantumespresso.workflows.protocols.utils import recursive_merge

    left_dict = {'a': {'b': 1, 'c': {'d': 2}}, 'g': 3}
    right_dict = {'a': {'c': {'d': 'D'}}, 'e': {'f': 'F'}}
    merged = recursive_merge(left_dict, right_dict)

    assert right_dict == {'a': {'c': {'d': 'D'}}, 'e': {'f': 'F'}}

    assert merged == {'a': {'b': 1, 'c': {'d': 'D'}}, 'e': {'f': 'F'}, 'g': 3}
