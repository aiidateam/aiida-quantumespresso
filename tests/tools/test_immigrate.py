# -*- coding: utf-8 -*-
# pylint: disable=unused-argument
"""Tests for immigrating `PwCalculation`s."""
from __future__ import absolute_import
import os

from aiida.common.folders import SandboxFolder
from aiida.orm import Dict, RemoteData

from aiida_quantumespresso.tools.pwinputparser import create_builder_from_file
from aiida_quantumespresso.tools.immigrate_calc import immigrate_existing


def test_create_builder(fixture_database, fixture_computer_localhost, generate_code_localhost, generate_upf_data,
                        generate_calc_job, fixture_sandbox_folder):
    """ this test uses the input file generated from tests.calculations.test_pw.test_pw_default"""
    entry_point_name = 'quantumespresso.pw'
    code = generate_code_localhost(entry_point_name, fixture_computer_localhost)

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

    in_foldername = os.path.join('tests', 'calculations', 'test_pw')
    in_folderpath = os.path.abspath(in_foldername)

    upf_foldername = os.path.join('tests', 'fixtures', 'pseudos')
    upf_folderpath = os.path.abspath(upf_foldername)
    si_upf = generate_upf_data('Si')
    si_upf.store()

    builder = create_builder_from_file(in_folderpath, 'test_pw_default.in', code, metadata, upf_folderpath)

    assert builder['code'] == code
    assert builder['metadata'] == metadata
    assert builder['pseudos']['Si'].id == si_upf.id
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

    generate_calc_job(fixture_sandbox_folder, entry_point_name, builder)


def test_immigrate_calc_pw(fixture_database, fixture_computer_localhost,
                           generate_code_localhost, generate_upf_data, data_regression):
    """ this test uses the input file generated from tests.calculations.test_pw.test_pw_default,
    and input file generated from tests.parsers.test_pw.test_pw_default
    """
    entry_point_name = 'quantumespresso.pw'
    code = generate_code_localhost(entry_point_name, fixture_computer_localhost)

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

    in_foldername = os.path.join('tests', 'calculations', 'test_pw')
    in_folderpath = os.path.abspath(in_foldername)

    upf_foldername = os.path.join('tests', 'fixtures', 'pseudos')
    upf_folderpath = os.path.abspath(upf_foldername)
    si_upf = generate_upf_data('Si')
    si_upf.store()

    out_foldername = os.path.join('tests', 'parsers', 'fixtures', 'pw', 'default')
    out_folderpath = os.path.abspath(out_foldername)

    builder = create_builder_from_file(in_folderpath, 'test_pw_default.in', code, metadata, upf_folderpath)
    builder.settings = Dict(dict={'NO_BANDS': True})
    fixture_computer_localhost.configure()

    with SandboxFolder() as folder:
        folder.insert_path(os.path.join(out_folderpath, 'aiida.out'))
        outfolder = folder.get_subfolder('out/aiida.save/', create=True)
        outfolder.insert_path(os.path.join(out_folderpath, 'data-file.xml'))

        remote = RemoteData(remote_path=folder.abspath, computer=fixture_computer_localhost)
        calc_node = immigrate_existing(builder, remote)

    attributes = calc_node.attributes
    attributes['remote_workdir'] = 'path/to/remote'
    attributes.pop('retrieve_singlefile_list', None)  # removed post v1.0.0b4
    attributes.pop('version', None)

    data_regression.check(attributes)

    outputs = calc_node.get_outgoing()
    assert set(outputs.all_link_labels()) == set(
        ['retrieved', 'remote_folder', 'output_trajectory',
         'output_parameters', 'output_kpoints'])

    assert calc_node.outputs.output_parameters.get_dict()['parser_warnings'] == []
