# -*- coding: utf-8 -*-
"""Tests for the ``XspectraCoreWorkChain.get_builder_from_protocol`` method."""
import io

from aiida.engine import ProcessBuilder
from aiida.orm import SinglefileData

from aiida_quantumespresso.workflows.xspectra.core import XspectraCoreWorkChain


def test_get_available_protocols():
    """Test ``XspectraCoreWorkChain.get_available_protocols``."""
    protocols = XspectraCoreWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['fast', 'moderate', 'precise']
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``XspectraCoreWorkChain.get_default_protocol``."""
    assert XspectraCoreWorkChain.get_default_protocol() == 'moderate'


def test_get_available_treatments():
    """Test ``XspectraCoreWorkChain.get_available_treatments``."""
    treatments = XspectraCoreWorkChain.get_available_treatments()
    assert sorted(treatments.keys()) == ['full', 'half', 'none', 'xch_fixed', 'xch_smear']
    assert all('description' in treatment for treatment in treatments.values())


def test_get_default_treatment():
    """Test ``XspectraCoreWorkChain.get_default_treatment``."""
    assert XspectraCoreWorkChain.get_default_treatment() == 'full'


def test_default(fixture_code, generate_structure, data_regression, serialize_builder):
    """Test ``XspectraCoreWorkChain.get_builder_from_protocol`` for the default protocol."""
    pw_code = fixture_code('quantumespresso.pw')
    xs_code = fixture_code('quantumespresso.xspectra')
    structure = generate_structure('silicon')
    core_wfc_data = SinglefileData(
        io.StringIO(
            '# number of core states 3 =  1 0;  2 0;'
            '\n6.51344e-05 6.615743462459999e-3'
            '\n6.59537e-05 6.698882211449999e-3'
        )
    )
    overrides = {'abs_atom_marker': 'Si'}
    builder = XspectraCoreWorkChain.get_builder_from_protocol(
        pw_code=pw_code, xs_code=xs_code, core_wfc_data=core_wfc_data, structure=structure, overrides=overrides
    )

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))
