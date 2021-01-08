# -*- coding: utf-8 -*-
"""Tests for the `create_magnetic_allotrope` calculation function."""
from aiida.orm import Float, List
from aiida.plugins import CalculationFactory
import pytest

create_magnetic_allotrope = CalculationFactory('quantumespresso.create_magnetic_allotrope')


@pytest.mark.usefixtures('aiida_profile')
def test_configuration_00(generate_structure_from_kinds):
    """Test `create_magnetic_allotrope` calculation function.

    Case: one kind but with equal magnetic moments.
    Expected result: no new kind names should be introduced.
    """
    kind_names = ['Fe', 'Fe']
    magnetic_moments = List(list=[0.2, 0.2])

    structure = generate_structure_from_kinds(kind_names)
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments).values()

    assert set(allotrope.get_kind_names()) == {'Fe'}
    assert allotrope_magnetic_moments.get_dict() == {'Fe': 0.2}


@pytest.mark.usefixtures('aiida_profile')
def test_configuration_01(generate_structure_from_kinds):
    """Test `create_magnetic_allotrope` calculation function.

    Case: two kinds all with equal magnetic moments.
    Expected result: no new kind names should be introduced.
    """
    kind_names = ['Fe', 'Fe', 'Ni', 'Ni', 'Ni']
    magnetic_moments = List(list=[0.2, 0.2, 0.5, 0.5, 0.5])

    structure = generate_structure_from_kinds(kind_names)
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments).values()

    assert set(allotrope.get_kind_names()) == {'Fe', 'Ni'}
    assert allotrope_magnetic_moments.get_dict() == {'Fe': 0.2, 'Ni': 0.5}


@pytest.mark.usefixtures('aiida_profile')
def test_configuration_02(generate_structure_from_kinds):
    """Test `create_magnetic_allotrope` calculation function.

    Case: only one kind but with unequal magnetic moments.
    Expected result: two new kinds introduced one for each magnetic moment.
    """
    kind_names = ['Fe', 'Fe']
    magnetic_moments = List(list=[0.2, 1.0])

    structure = generate_structure_from_kinds(kind_names)
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments).values()

    assert set(allotrope.get_kind_names()) == {'Fe0', 'Fe1'}
    assert allotrope_magnetic_moments.get_dict() == {'Fe0': 1.0, 'Fe1': 0.2}


@pytest.mark.usefixtures('aiida_profile')
def test_configuration_03(generate_structure_from_kinds):
    """Test `create_magnetic_allotrope` calculation function.

    Case: only one kind but with three types of magnetic moments that are not grouped together.
    Expected result: two new kinds introduced one for each magnetic moment.
    """
    kind_names = ['Fe', 'Fe', 'Fe', 'Fe']
    magnetic_moments = List(list=[0.2, 0.8, 1.5, 0.8])

    structure = generate_structure_from_kinds(kind_names)
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments).values()

    assert set(allotrope.get_kind_names()) == {'Fe0', 'Fe1', 'Fe2'}
    assert allotrope_magnetic_moments.get_dict() == {'Fe0': 1.5, 'Fe1': 0.8, 'Fe2': 0.2}


@pytest.mark.usefixtures('aiida_profile')
def test_configuration_04(generate_structure_from_kinds):
    """Test `create_magnetic_allotrope` calculation function.

    Case: only one kind but with four different values of magnetic moments but middle two are within tolerance.
    Expected result: two new kinds introduced one for each magnetic moment.
    """
    kind_names = ['Fe', 'Fe', 'Fe', 'Fe']
    magnetic_moments = List(list=[0.0, 0.50, 0.45, 0.40])

    structure = generate_structure_from_kinds(kind_names)

    # Default tolerances: just two different kinds
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments).values()
    assert set(allotrope.get_kind_names()) == {'Fe0', 'Fe1'}
    assert [site.kind_name for site in allotrope.sites] == ['Fe0', 'Fe1', 'Fe1', 'Fe1']
    assert allotrope_magnetic_moments.get_dict() == {'Fe0': 0.0, 'Fe1': 0.5}

    # Lower atol to 0.05: 0.5 & 0.45 now one kind, 0.4 new kind -> three different kinds
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments,
                                                                      atol=Float(0.05)).values()
    assert set(allotrope.get_kind_names()) == {'Fe0', 'Fe1', 'Fe2'}
    assert [site.kind_name for site in allotrope.sites] == ['Fe0', 'Fe1', 'Fe1', 'Fe2']
    assert allotrope_magnetic_moments.get_dict() == {'Fe0': 0.0, 'Fe1': 0.5, 'Fe2': 0.4}

    # Increase atol to 0.1, again only two different kinds
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments,
                                                                      atol=Float(0.1)).values()
    assert set(allotrope.get_kind_names()) == {'Fe0', 'Fe1'}
    assert [site.kind_name for site in allotrope.sites] == ['Fe0', 'Fe1', 'Fe1', 'Fe1']
    assert allotrope_magnetic_moments.get_dict() == {'Fe0': 0.0, 'Fe1': 0.5}

    # Really strict tolerance or atol = 0.01: All sites get different kinds
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments,
                                                                      atol=Float(1E-2)).values()
    assert set(allotrope.get_kind_names()) == {'Fe0', 'Fe1', 'Fe2', 'Fe3'}
    assert [site.kind_name for site in allotrope.sites] == ['Fe0', 'Fe1', 'Fe2', 'Fe3']
    assert allotrope_magnetic_moments.get_dict() == {'Fe0': 0.0, 'Fe1': 0.5, 'Fe2': 0.45, 'Fe3': 0.4}


@pytest.mark.usefixtures('aiida_profile')
def test_configuration_05(generate_structure_from_kinds):
    """Test `create_magnetic_allotrope` calculation function.

    Case: One kind, only negative magnetic moments with one close to zero
    Expected result: Depends on tolerance, see below
    """
    kind_names = ['Fe', 'Fe', 'Fe', 'Fe']
    magnetic_moments = List(list=[-0.5, -0.6, -1.5, -0.01])

    structure = generate_structure_from_kinds(kind_names)

    # Default tolerance values, one zero site and two magnetic
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments).values()
    assert set(allotrope.get_kind_names()) == {'Fe0', 'Fe1', 'Fe2'}
    assert [site.kind_name for site in allotrope.sites] == ['Fe0', 'Fe1', 'Fe2', 'Fe2']
    assert allotrope_magnetic_moments.get_dict() == {'Fe0': 0.0, 'Fe1': -1.5, 'Fe2': -0.6}

    # Strict absolute tolerance, one zero site and three magnetic
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments,
                                                                      atol=Float(0.05)).values()
    assert set(allotrope.get_kind_names()) == {'Fe0', 'Fe1', 'Fe2', 'Fe3'}
    assert [site.kind_name for site in allotrope.sites] == ['Fe0', 'Fe1', 'Fe2', 'Fe3']
    assert allotrope_magnetic_moments.get_dict() == {'Fe0': 0.0, 'Fe1': -1.5, 'Fe2': -0.6, 'Fe3': -0.5}

    # Strict absolute and zero tolerance, four magnetic sites
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(
        structure, magnetic_moments, atol=Float(0.05), ztol=Float(1E-3)
    ).values()
    assert set(allotrope.get_kind_names()) == {'Fe0', 'Fe1', 'Fe2', 'Fe3'}
    assert [site.kind_name for site in allotrope.sites] == ['Fe0', 'Fe1', 'Fe2', 'Fe3']
    assert allotrope_magnetic_moments.get_dict() == {'Fe0': -1.5, 'Fe1': -0.6, 'Fe2': -0.5, 'Fe3': -0.01}


@pytest.mark.usefixtures('aiida_profile')
def test_configuration_06(generate_structure_from_kinds):
    """Test `create_magnetic_allotrope` calculation function.

    Case: Two kinds, magnetic moments with different signs for the first (Fe)
    Expected result: Depends on tolerance, see below
    """
    kind_names = ['Fe', 'Fe', 'Fe', 'Fe', 'Ni', 'Ni']
    magnetic_moments = List(list=[-0.1, 0.1, -0.2, 0.01, 0.2, 0.25])

    structure = generate_structure_from_kinds(kind_names)

    # Default tolerance values, one zero and two magnetic sites for Fe, one magnetic site for Ni
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments).values()
    assert set(allotrope.get_kind_names()) == {'Fe0', 'Fe1', 'Fe2', 'Ni'}
    assert allotrope_magnetic_moments.get_dict() == {'Fe0': 0.0, 'Fe1': 0.1, 'Fe2': -0.2, 'Ni': 0.25}

    # Very strict absolute tolerance, all different sites
    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments,
                                                                      atol=Float(0.02)).values()
    assert set(allotrope.get_kind_names()) == {'Fe0', 'Fe1', 'Fe2', 'Fe3', 'Ni0', 'Ni1'}
    assert allotrope_magnetic_moments.get_dict() == {
        'Fe0': 0.0,
        'Fe1': 0.1,
        'Fe2': -0.2,
        'Fe3': -0.1,
        'Ni0': 0.25,
        'Ni1': 0.2
    }


@pytest.mark.usefixtures('aiida_profile')
def test_configuration_07(generate_structure_from_kinds):
    """Test `create_magnetic_allotrope` calculation function.

    Case: Two different symbols but the same magnetic moment.
    Expected result: One kind with name equal to the element symbol
    """
    kind_names = ['Fe0', 'Fe1']
    magnetic_moments = List(list=[0.1, 0.1])

    structure = generate_structure_from_kinds(kind_names)

    allotrope, allotrope_magnetic_moments = create_magnetic_allotrope(structure, magnetic_moments).values()
    assert set(allotrope.get_kind_names()) == {'Fe'}
    assert allotrope_magnetic_moments.get_dict() == {'Fe': 0.1}
