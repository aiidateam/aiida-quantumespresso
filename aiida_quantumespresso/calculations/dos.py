# -*- coding: utf-8 -*-
from __future__ import absolute_import

from aiida import orm
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class DosCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the dos.x code of Quantum ESPRESSO."""

    _DOS_FILENAME = 'aiida.dos'
    _default_namelists = ['DOS']
    _blocked_keywords = [
        ('DOS', 'fildos', _DOS_FILENAME),
        ('DOS', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),
        ('DOS', 'prefix', NamelistsCalculation._PREFIX),
    ]
    _internal_retrieve_list = [_DOS_FILENAME]
    _default_parser = 'quantumespresso.dos'

    @classmethod
    def define(cls, spec):
        super(DosCalculation, cls).define(spec)
        spec.input('parent_folder', valid_type=(orm.RemoteData, orm.FolderData), required=True)
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('output_dos', valid_type=orm.XyData)
        spec.default_output_node = 'output_parameters'
        spec.exit_code(
            100, 'ERROR_NO_RETRIEVED_FOLDER', message='The retrieved folder data node could not be accessed.')
        spec.exit_code(
            110, 'ERROR_READING_OUTPUT_FILE', message='The output file could not be read from the retrieved folder.')
        spec.exit_code(
            111, 'ERROR_READING_DOS_FILE', message='The dos file could not be read from the retrieved folder.')
        spec.exit_code(
            130, 'ERROR_JOB_NOT_DONE', message='The computation did not finish properly ("JOB DONE" not found).')
