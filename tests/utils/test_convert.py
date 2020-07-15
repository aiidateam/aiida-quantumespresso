# -*- coding: utf-8 -*-
"""Tests for :py:mod:`~aiida_quantumespresso.utils.convert`."""
import unittest

from aiida_quantumespresso.utils.convert import convert_input_to_namelist_entry


class TestUtilsConvertInputToNamelistEntry(unittest.TestCase):
    """Tests for :py:func:`~aiida_quantumespresso.utils.convert.convert_input_to_namelist_entry`."""

    @classmethod
    def setUpClass(cls):
        """Define a test mapping."""
        super().setUpClass()
        cls.mapping = {
            'Co': 1,
            'O': 3,
        }

    def validate_converted_output(self, parameters, expected):
        """Validate recursively the two dictionaries."""
        for key, value in parameters.items():
            converted = convert_input_to_namelist_entry(key, value, self.mapping)
            lines = [line for line in converted.split('\n') if line]

            for index, line in enumerate(lines):
                self.assertIn(expected[key][index], line)

    def test_simple_value(self):
        """When the value is a single value, it should be formatted as 'key = value'."""
        parameters = {'ecutwfc': 100}
        expected = {'ecutwfc': ['ecutwfc = 100']}
        self.validate_converted_output(parameters, expected)

    def test_single_list(self):
        """When the value is a single list, each entry should be formatted as 'key(i) = value[i]'."""
        parameters = {'efield': [4, 5, 6]}
        expected = {
            'efield': [
                'efield(1) = 4',
                'efield(2) = 5',
                'efield(3) = 6',
            ]
        }
        self.validate_converted_output(parameters, expected)

    def test_double_list(self):
        """For a double nested list, each sub list will be formatted as 'key(values[i],..,values[n-1]) = values[n]'."""
        parameters = {
            'starting_ns_eigenvalues': [
                [1, 1, 2, 10],
                [2, 2, 3, 10],
            ]
        }
        expected = {
            'starting_ns_eigenvalues': [
                'starting_ns_eigenvalues(1,1,2) = 10',
                'starting_ns_eigenvalues(2,2,3) = 10',
            ]
        }
        self.validate_converted_output(parameters, expected)

    def test_double_list_kind_mapping(self):
        """If any value in the sub lists matches a kind name it should be mapped using the mapping."""
        parameters = {
            'hubbard_j': [
                [1, 'Co', 10],
                [2, 'O', 10],
            ]
        }
        expected = {
            'hubbard_j': [
                'hubbard_j(1,1) = 10',
                'hubbard_j(2,3) = 10',
            ]
        }
        self.validate_converted_output(parameters, expected)

    def test_dictionary(self):
        """For a dictionary each key should be mapped using the mapping: 'key(mapping[k]) = value'."""
        parameters = {
            'hubbard_u': {
                'Co': 1,
                'O': 2,
            },
        }
        expected = {
            'hubbard_u': [
                'hubbard_u(1) = 1',
                'hubbard_u(3) = 2',
            ]
        }
        self.validate_converted_output(parameters, expected)

    def test_double_list_invalid_type(self):
        """Values in the double nested lists should be integers or strings."""
        parameters = {
            'starting_ns_eigenvalues': [
                [1, (), 2, 10],
            ]
        }

        for key, value in parameters.items():
            with self.assertRaises(ValueError):
                convert_input_to_namelist_entry(key, value, self.mapping)

    def test_double_list_invalid_kind(self):
        """If a value in the nested list is a string it should be present as a key in the mapping dictionary."""
        parameters = {
            'starting_ns_eigenvalues': [
                [1, 'Ni', 2, 10],
            ]
        }

        for key, value in parameters.items():
            with self.assertRaises(ValueError):
                convert_input_to_namelist_entry(key, value, self.mapping)

    def test_double_list_no_mapping(self):
        """If the double list contains valid kinds, a mapping has to be provided."""
        parameters = {
            'starting_ns_eigenvalues': [
                [1, 'Co', 2, 10],
            ]
        }

        for key, value in parameters.items():
            with self.assertRaises(ValueError):
                convert_input_to_namelist_entry(key, value, None)
