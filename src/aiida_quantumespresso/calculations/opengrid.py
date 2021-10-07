# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the open_grid.x code of Quantum ESPRESSO."""
from aiida.orm import RemoteData, FolderData, KpointsData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class OpengridCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the open_grid.x code of Quantum ESPRESSO."""

    _default_namelists = ['INPUTPP']
    _blocked_keywords = [
        ('INPUTPP', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),
        ('INPUTPP', 'prefix', NamelistsCalculation._PREFIX),
        ('INPUTPP', 'overwrite_prefix', True),
    ]
    _default_parser = 'quantumespresso.opengrid'

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('parent_folder', valid_type=(RemoteData, FolderData),
            help='The output folder of a completed `PwCalculation` on an irreducible Brillouin zone')
        # Outputs the dimensions of the unfolded kmesh
        spec.output('kpoints_mesh', valid_type=KpointsData)
        # Outputs an explicit list of kpoints of the unfolded kmesh
        spec.output('kpoints', valid_type=KpointsData)

        spec.exit_code(300, 'ERROR_NO_RETRIEVED_FOLDER',
            message='The retrieved folder data node could not be accessed.')
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(311, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete probably because the calculation got interrupted.')
        spec.exit_code(312, 'ERROR_INCOMPATIBLE_FFT_GRID',
            message='Found rotation or fractional translation not compatible with FFT grid.')
        spec.exit_code(340, 'ERROR_GENERIC_QE_ERROR',
            message='Encountered a generic error message.')
        spec.exit_code(350, 'ERROR_OUTPUT_KPOINTS_MISMATCH',
            message='Mismatch between kmesh dimensions and number of kpoints.')
