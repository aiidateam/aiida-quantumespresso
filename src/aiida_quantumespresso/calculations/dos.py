# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the dos.x code of Quantum ESPRESSO."""
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
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('parent_folder', valid_type=(orm.RemoteData, orm.FolderData), required=True)
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('output_dos', valid_type=orm.XyData)
        spec.default_output_node = 'output_parameters'
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete probably because the calculation got interrupted.')
        spec.exit_code(330, 'ERROR_READING_DOS_FILE',
            message='The dos file could not be read from the retrieved folder.')
        # yapf: enable
