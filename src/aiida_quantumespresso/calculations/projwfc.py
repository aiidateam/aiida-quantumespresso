"""`CalcJob` implementation for the projwfc.x code of Quantum ESPRESSO."""

from pathlib import Path

from aiida.orm import Dict, FolderData, RemoteData, XyData

from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class ProjwfcCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the projwfc.x code of Quantum ESPRESSO.

    Projwfc.x code of the Quantum ESPRESSO distribution, handles the the computation of projections of bloch
    wavefunctions onto atomic orbitals.

    <Psi(n,k) | Y(theta,phi)R(r) >. For more information, refer to http://www.quantum-espresso.org/
    """

    _default_namelists = ['PROJWFC']
    _blocked_keywords = [
        ('PROJWFC', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),  # noqa: SLF001
        ('PROJWFC', 'prefix', NamelistsCalculation._PREFIX),  # noqa: SLF001
        ('PROJWFC', 'lsym', True),
        ('PROJWFC', 'lwrite_overlaps', False),
        ('PROJWFC', 'lbinary_data', False),
        ('PROJWFC', 'kresolveddos', False),
    ]
    _default_parser = 'quantumespresso.projwfc'

    xml_path = Path(NamelistsCalculation._default_parent_output_folder).joinpath(  # noqa: SLF001
        f'{NamelistsCalculation._PREFIX}.save',  # noqa: SLF001
        'data-file-schema.xml',
    )

    # The XML file is added to the temporary retrieve list since it is required for parsing, but already in the
    # repository of a an ancestor calculation.
    _retrieve_temporary_list = [
        NamelistsCalculation._PREFIX + '.pdos*',  # noqa: SLF001
        NamelistsCalculation._PREFIX + '.ldos_boxes',  # noqa: SLF001
        xml_path.as_posix(),
    ]

    @classmethod
    def define(cls, spec):
        """Define the process specification."""

        from aiida.orm import BandsData, ProjectionData

        super().define(spec)
        spec.input('parent_folder', valid_type=(RemoteData, FolderData), help='The output folder of a pw.x calculation')
        spec.output('output_parameters', valid_type=Dict)
        spec.output('Dos', valid_type=XyData, help='Total DOS')
        # if tdosinboxes: an ldos_boxes file is generated
        spec.output('Ldos', valid_type=XyData, required=False, help='LDOS for each box in same XyData node')
        # if not tdosinboxes: no pdos_tot file is generated
        spec.output(
            'Pdos',
            valid_type=XyData,
            required=False,
            help='Total Projected DOS (on all orbitals or on all boxes if LDOS)',
        )
        # if spin: Dos and Pdos have a second y-array for the spin down
        spec.output('projections_up', valid_type=ProjectionData, required=False)
        spec.output('projections_down', valid_type=ProjectionData, required=False)
        spec.output('bands_up', valid_type=BandsData, required=False)
        spec.output('bands_down', valid_type=BandsData, required=False)
        # if non-spin
        spec.output('projections', valid_type=ProjectionData, required=False)
        spec.output('bands', valid_type=BandsData, required=False)

        spec.default_output_node = 'output_parameters'

        spec.exit_code(
            301, 'ERROR_NO_RETRIEVED_TEMPORARY_FOLDER', message='The retrieved temporary folder could not be accessed.'
        )
        spec.exit_code(
            302,
            'ERROR_OUTPUT_STDOUT_MISSING',
            message='The retrieved folder did not contain the required stdout output file.',
        )
        spec.exit_code(
            303, 'ERROR_OUTPUT_XML_MISSING', message='The retrieved folder did not contain the required XML file.'
        )
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ', message='The stdout output file could not be read.')
        spec.exit_code(311, 'ERROR_OUTPUT_STDOUT_PARSE', message='The stdout output file could not be parsed.')
        spec.exit_code(
            312,
            'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete probably because the calculation got interrupted.',
        )
        spec.exit_code(320, 'ERROR_OUTPUT_XML_READ', message='The XML output file could not be read.')
        spec.exit_code(321, 'ERROR_OUTPUT_XML_PARSE', message='The XML output file could not be parsed.')
        spec.exit_code(322, 'ERROR_OUTPUT_XML_FORMAT', message='The XML output file has an unsupported format.')
        spec.exit_code(
            330, 'ERROR_READING_PDOSTOT_FILE', message='The pdos_tot file could not be read from the retrieved folder.'
        )
        spec.exit_code(
            331,
            'ERROR_READING_LDOSBOXES_FILE',
            message='The ldos_boxes file could not be read from the retrieved folder.',
        )
        spec.exit_code(
            332, 'ERROR_MISSING_PDOSTOT_FILE', message='The pdos_tot file is missing from the retrieved folder.'
        )
        spec.exit_code(
            333, 'ERROR_MISSING_LDOSBOXES_FILE', message='The ldos_boxes file is missing from the retrieved folder.'
        )
        spec.exit_code(
            340, 'ERROR_PARSING_PROJECTIONS', message='An exception was raised parsing bands and projections.'
        )
        spec.exit_code(
            350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION', message='The parser raised an unexpected exception: {exception}'
        )
