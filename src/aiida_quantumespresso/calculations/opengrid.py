# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the open_grid.x code of Quantum ESPRESSO."""
from aiida.orm import RemoteData, FolderData, Dict, KpointsData, StructureData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class OpengridCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the open_grid.x code of Quantum ESPRESSO."""

    _default_namelists = ['INPUTPP']
    _blocked_keywords = [
        ('INPUTPP', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),
        ('INPUTPP', 'prefix', NamelistsCalculation._PREFIX),
    ]
    _default_parser = 'quantumespresso.opengrid'

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('parent_folder', valid_type=(RemoteData, FolderData),
                   help='The output folder of a completed `PwCalculation`')
        # ProjwfcParser of ProjwfcCalculation requires a previous Calculation
        # has an input or output structure to parse atomic orbitals.
        # Here defines an optional input structure so ProjwfcCalculation
        # could restart from an OpengridCalculation.
        # We cannot output the structure of the parent calculation, otherwise
        # the following exception occurs for the StructureData node
        #   ValueError: node already has an incoming LinkType.CREATE link
        spec.input('structure', valid_type=StructureData, required=False,
                    help='Could be used by a subsequent projwfc.x calculation')
        # ProjwfcParser also requires input or output kpoints,
        # output_parameters which must has number_of_spin_components
        spec.output('output_parameters', valid_type=Dict,
                    help='Containing number_of_spin_components, could be used by a subsequent projwfc.x calculation')
        # kpoints_mesh outputs the dimensions of unfolded kmesh
        spec.output('kpoints_mesh', valid_type=KpointsData,
                    help='The dimensions of unfolded kpoints mesh')
        # kpoints outputs the explicit points of unfolded kmesh
        spec.output('kpoints', valid_type=KpointsData,
                    help='An explicit list of unfolded kpoints')
        spec.default_output_node = 'output_parameters'
        spec.exit_code(300, 'ERROR_NO_RETRIEVED_FOLDER',
            message='The retrieved folder data node could not be accessed.')
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete probably because the calculation got interrupted.')
        spec.exit_code(340, 'ERROR_GENERIC_QE_ERROR',
            message='Encountered a generic error message')
        spec.exit_code(350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION',
            message='An error happened while parsing the output file')
