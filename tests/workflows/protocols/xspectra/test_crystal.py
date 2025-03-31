# -*- coding: utf-8 -*-
"""Tests for the ``XspectraCrystalWorkChain.get_builder_from_protocol`` method."""
import io

from aiida.engine import ProcessBuilder
from aiida.orm import SinglefileData

from aiida_quantumespresso.workflows.xspectra.crystal import XspectraCrystalWorkChain


def test_get_available_protocols():
    """Test ``XspectraCrystalWorkChain.get_available_protocols``."""
    protocols = XspectraCrystalWorkChain.get_available_protocols()
    assert sorted(protocols.keys()) == ['balanced', 'fast', 'stringent']
    assert all('description' in protocol for protocol in protocols.values())


def test_get_default_protocol():
    """Test ``XspectraCrystalWorkChain.get_default_protocol``."""
    assert XspectraCrystalWorkChain.get_default_protocol() == 'balanced'


def test_default(fixture_code, generate_structure, data_regression, serialize_builder, generate_upf_data):
    """Test ``XspectraCrystalWorkChain.get_builder_from_protocol`` for the default protocol."""
    pw_code = fixture_code('quantumespresso.pw')
    xs_code = fixture_code('quantumespresso.xspectra')
    structure = generate_structure('silicon')
    pseudos = {'Si': {'core_hole': generate_upf_data('Si'), 'gipaw': generate_upf_data('Si')}}
    core_wfc_data = {
        'Si':
        SinglefileData(
            io.StringIO(
                '# number of core states 3 =  1 0;  2 0;'
                '\n6.51344e-05 6.615743462459999e-3'
                '\n6.59537e-05 6.698882211449999e-3'
            )
        )
    }
    overrides = {'abs_atom_marker': 'Si'}
    builder = XspectraCrystalWorkChain.get_builder_from_protocol(
        pw_code=pw_code,
        xs_code=xs_code,
        core_wfc_data=core_wfc_data,
        structure=structure,
        overrides=overrides,
        pseudos=pseudos
    )

    assert isinstance(builder, ProcessBuilder)
    data_regression.check(serialize_builder(builder))
