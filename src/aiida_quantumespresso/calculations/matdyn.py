# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the matdyn.x code of Quantum ESPRESSO."""
from aiida import orm

from aiida_quantumespresso.calculations.namelists import NamelistsCalculation
from aiida_quantumespresso.data.force_constants import ForceConstantsData


class MatdynCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the matdyn.x code of Quantum ESPRESSO."""

    _PHONON_FREQUENCIES_NAME = 'phonon_frequencies.dat'
    _PHONON_MODES_NAME = 'phonon_displacements.dat'
    _PHONON_DOS_NAME = 'phonon_dos.dat'

    _default_namelists = ['INPUT']
    _blocked_keywords = [
        ('INPUT', 'flfrq', _PHONON_FREQUENCIES_NAME),  # output frequencies
        ('INPUT', 'flvec', _PHONON_MODES_NAME),  # output displacements
        ('INPUT', 'fldos', _PHONON_DOS_NAME),  # output density of states
        ('INPUT', 'q_in_cryst_coord', True),  # kpoints always in crystal coordinates
    ]

    _internal_retrieve_list = [_PHONON_FREQUENCIES_NAME, _PHONON_DOS_NAME]
    _default_parser = 'quantumespresso.matdyn'

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('force_constants', valid_type=ForceConstantsData, required=True)
        spec.input('kpoints', valid_type=orm.KpointsData, help='Kpoints on which to calculate the phonon frequencies.')
        spec.inputs.validator = cls._validate_inputs

        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('output_phonon_bands', valid_type=orm.BandsData)
        spec.default_output_node = 'output_parameters'

        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(312, 'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete probably because the calculation got interrupted.')
        spec.exit_code(330, 'ERROR_OUTPUT_FREQUENCIES',
            message='The output frequencies file could not be read from the retrieved folder.')
        spec.exit_code(410, 'ERROR_OUTPUT_KPOINTS_MISSING',
            message='Number of kpoints not found in the output data')
        spec.exit_code(411, 'ERROR_OUTPUT_KPOINTS_INCOMMENSURATE',
            message='Number of kpoints in the inputs is not commensurate with those in the output')
        # yapf: enable

    @staticmethod
    def _validate_inputs(value, _):
        """Validate the top level namespace."""
        if 'parameters' in value:
            parameters = value['parameters'].get_dict()
        else:
            parameters = {}

        if parameters.get('INPUT', {}).get('flfrc', None) is not None:
            return '`INPUT.flfrc` is set automatically from the `force_constants` input.'

    def generate_input_file(self, parameters):
        """Generate namelist input_file content given a dict of parameters.

        :param parameters: 'dict' containing the fortran namelists and parameters to be used.
          e.g.: {'CONTROL':{'calculation':'scf'}, 'SYSTEM':{'ecutwfc':30}}
        :return: 'str' containing the input_file content a plain text.
        """
        kpoints = self.inputs.kpoints
        parameters.setdefault('INPUT', {})['flfrc'] = self.inputs.force_constants.filename
        file_content = super().generate_input_file(parameters)

        try:
            kpoints_list = kpoints.get_kpoints()
        except AttributeError:
            kpoints_list = kpoints.get_kpoints_mesh(print_list=True)

        kpoints_string = [f'{len(kpoints_list)}']
        for kpoint in kpoints_list:
            kpoints_string.append('{:18.10f} {:18.10f} {:18.10f}'.format(*kpoint))  # pylint: disable=consider-using-f-string

        file_content += '\n'.join(kpoints_string) + '\n'

        return file_content

    def prepare_for_submission(self, folder):
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :class:`~aiida.common.datastructures.CalcInfo` instance.
        """
        calcinfo = super().prepare_for_submission(folder)

        force_constants = self.inputs.force_constants
        calcinfo.local_copy_list.append((force_constants.uuid, force_constants.filename, force_constants.filename))

        return calcinfo
