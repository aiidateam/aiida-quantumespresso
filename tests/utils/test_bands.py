# -*- coding: utf-8 -*-
"""Unit tests for the :py:mod:`~aiida_quantumespresso.utils.bands` module."""
from __future__ import absolute_import

import numpy
import pytest

from aiida_quantumespresso.utils.bands import get_highest_occupied_band


class TestGetHighestOccupiedBand(object):
    """Tests for :py:func:`~aiida_quantumespresso.utils.bands.get_highest_occupied_band`."""

    def test_valid_node(self, fixture_database):
        """Test that the correct exceptions are thrown for incompatible nodes."""
        from aiida.orm import ArrayData, BandsData

        # Invalid node type
        node = ArrayData().store()
        with pytest.raises(ValueError):
            get_highest_occupied_band(node)

        # The `occupations` array is missing
        node = BandsData()
        node.set_array('not_occupations', numpy.array([]))
        node.store()
        with pytest.raises(ValueError):
            get_highest_occupied_band(node)

        # The `occupations` array has incorrect shape
        node = BandsData()
        node.set_array('occupations', numpy.array([1., 1.]))
        node.store()
        with pytest.raises(ValueError):
            get_highest_occupied_band(node)

    def test_spin_unpolarized(self, fixture_database):
        """Test the function for a non spin-polarized calculation meaning there will be a single spin channel."""
        from aiida.orm import BandsData

        occupations = numpy.array([
            [2., 2., 2., 2., 0.],
            [2., 2., 2., 2., 0.],
            [2., 2., 2., 2., 0.],
            [2., 2., 2., 2., 0.],
        ])

        bands = BandsData()
        bands.set_array('occupations', occupations)
        bands.store()
        homo = get_highest_occupied_band(bands)
        assert homo == 4

    def test_spin_polarized(self, fixture_database):
        """Test the function for a spin-polarized calculation meaning there will be two spin channels."""
        from aiida.orm import BandsData

        occupations = numpy.array([
            [
                [2., 2., 2., 2., 0.],
                [2., 2., 2., 2., 0.],
            ],
            [
                [2., 2., 2., 2., 0.],
                [2., 2., 2., 2., 0.],
            ]
        ])

        bands = BandsData()
        bands.set_array('occupations', occupations)
        bands.store()
        homo = get_highest_occupied_band(bands)
        assert homo == 4
