# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the projwfc.x code of Quantum ESPRESSO."""
from aiida.orm import RemoteData, FolderData, Dict, XyData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class ProjwfcCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the projwfc.x code of Quantum ESPRESSO.

    Projwfc.x code of the Quantum ESPRESSO distribution, handles the the computation of projections of bloch
    wavefunctions onto atomic orbitals.

    <Psi(n,k) | Y(theta,phi)R(r) >. For more information, refer to http://www.quantum-espresso.org/
    """

    _default_namelists = ['PROJWFC']
    _blocked_keywords = [
        ('PROJWFC', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),
        ('PROJWFC', 'prefix', NamelistsCalculation._PREFIX),
        ('PROJWFC', 'lsym', True),
        ('PROJWFC', 'lwrite_overlaps', False),
        ('PROJWFC', 'lbinary_data', False),
        ('PROJWFC', 'kresolveddos', False),
        ('PROJWFC', 'tdosinboxes', False),
        ('PROJWFC', 'plotboxes', False),
    ]
    _default_parser = 'quantumespresso.projwfc'
    _internal_retrieve_list = [NamelistsCalculation._PREFIX + '.pdos*']

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        from aiida.orm import ProjectionData, BandsData
        super().define(spec)
        spec.input('parent_folder', valid_type=(RemoteData, FolderData), help='The output folder of a pw.x calculation')
        spec.output('output_parameters', valid_type=Dict)
        spec.output('Dos', valid_type=XyData)
        # if spin
        spec.output('projections_up', valid_type=ProjectionData, required=False)
        spec.output('projections_down', valid_type=ProjectionData, required=False)
        spec.output('bands_up', valid_type=BandsData, required=False)
        spec.output('bands_down', valid_type=BandsData, required=False)
        # if non-spin
        spec.output('projections', valid_type=ProjectionData, required=False)
        spec.output('bands', valid_type=BandsData, required=False)
        spec.default_output_node = 'output_parameters'
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete probably because the calculation got interrupted.')
        spec.exit_code(330, 'ERROR_READING_PDOSTOT_FILE',
            message='The pdos_tot file could not be read from the retrieved folder.')
        spec.exit_code(340, 'ERROR_PARSING_PROJECTIONS',
            message='An exception was raised parsing bands and projections.')
        # yapf: enable
