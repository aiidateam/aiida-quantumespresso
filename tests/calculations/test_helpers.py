# -*- coding: utf-8 -*-
"""Tests for the calculation input helper utilities."""
import pytest

from aiida_quantumespresso.calculations.helpers import pw_input_helper, QEInputValidationError


def test_pw_helper_multidimensional(generate_structure):
    """Test the helper for parameters containing a multidimensional parameter."""
    structure = generate_structure()
    parameters = {
        'CONTROL': {
            'calculation': 'scf'
        },
        'SYSTEM': {
            'ecutwfc': 30,
            'starting_ns_eigenvalue': [[1, 2, 'Si', 4]],
            'hubbard_j': [[1, 'Si', 15.7]]
        }
    }

    result = pw_input_helper(parameters, structure, version='6.4')

    assert result['CONTROL'] == parameters['CONTROL']
    assert result['SYSTEM'] == parameters['SYSTEM']

    with pytest.raises(QEInputValidationError):
        parameters['SYSTEM']['hubbard_j'] = [[1, 2, 15.7]]  # Second element is not a kind name
        pw_input_helper(parameters, structure, version='6.4')

    with pytest.raises(QEInputValidationError):
        parameters['SYSTEM']['hubbard_j'] = [[1, 'Ge', 15.7]]  # Second element is a non-existing structure kind name
        pw_input_helper(parameters, structure, version='6.4')
