# -*- coding: utf-8 -*-
"""Tests for immigrating `PwCalculation`s."""
import os

from aiida_quantumespresso.tools.pwinputparser import create_builder_from_file


def test_create_builder(
    aiida_profile, fixture_sandbox, fixture_code, generate_upf_data, generate_calc_job, filepath_tests
):
    """Test the `create_builder_from_file` method that parses an existing `pw.x` folder into a process builder.

    The input file used is the one generated for `tests.calculations.test_pw.test_pw_default`.
    """
    entry_point_name = 'quantumespresso.pw'
    code = fixture_code(entry_point_name)

    metadata = {
        'options': {
            'resources': {
                'num_machines': 1,
                'num_mpiprocs_per_machine': 32,
            },
            'max_memory_kb': 1000,
            'max_wallclock_seconds': 60 * 60 * 12,
            'withmpi': True,
        }
    }

    in_folderpath = os.path.join(filepath_tests, 'calculations', 'test_pw')
    upf_folderpath = os.path.join(filepath_tests, 'fixtures', 'pseudos')

    si_upf = generate_upf_data('Si')
    si_upf.store()

    builder = create_builder_from_file(
        in_folderpath, 'test_pw_default.in', code, metadata, upf_folderpath, use_first=True
    )

    assert builder['code'] == code
    assert builder['metadata'] == metadata
    pseudo_hash = si_upf.get_hash()
    assert pseudo_hash is not None
    assert builder['pseudos']['Si'].get_hash() == pseudo_hash
    assert builder['parameters'].get_dict() == {
        'CONTROL': {
            'calculation': 'scf',
            'verbosity': 'high'
        },
        'SYSTEM': {
            'ecutrho': 240.0,
            'ecutwfc': 30.0
        }
    }
    assert 'kpoints' in builder
    assert 'structure' in builder

    generate_calc_job(fixture_sandbox, entry_point_name, builder)
