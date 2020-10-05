# -*- coding: utf-8 -*-
"""Unit tests for the :py:mod:`~aiida_quantumespresso.utils.bands` module."""
import numpy
import pytest

from aiida_quantumespresso.utils.bands import get_highest_occupied_band


class TestGetHighestOccupiedBand:
    """Tests for :py:func:`~aiida_quantumespresso.utils.bands.get_highest_occupied_band`."""

    @staticmethod
    def test_valid_node():
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

    @staticmethod
    def test_threshold():
        """Test the `threshold` parameter."""
        from aiida.orm import BandsData

        threshold = 0.002

        bands = BandsData()
        bands.set_array('occupations', numpy.array([[2., 2., 2., 2., 0.001, 0.0015]]))
        bands.store()

        # All bands above the LUMO (occupation of 0.001) are below `2 * threshold`
        homo = get_highest_occupied_band(bands, threshold=threshold)
        assert homo == 4

        bands = BandsData()
        bands.set_array('occupations', numpy.array([[2., 2., 2., 2., 0.001, 0.003]]))
        bands.store()

        # A band above the LUMO (occupation of 0.001) has an occupation above `2 * threshold`
        with pytest.raises(ValueError):
            get_highest_occupied_band(bands, threshold=threshold)

    @staticmethod
    def test_spin_unpolarized():
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

    @staticmethod
    def test_spin_polarized():
        """Test the function for a spin-polarized calculation meaning there will be two spin channels."""
        from aiida.orm import BandsData

        occupations = numpy.array([[
            [2., 2., 2., 2., 0.],
            [2., 2., 2., 2., 0.],
        ], [
            [2., 2., 2., 2., 0.],
            [2., 2., 2., 2., 0.],
        ]])

        bands = BandsData()
        bands.set_array('occupations', occupations)
        bands.store()
        homo = get_highest_occupied_band(bands)
        assert homo == 4
