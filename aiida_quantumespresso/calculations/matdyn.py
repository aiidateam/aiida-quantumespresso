# -*- coding: utf-8 -*-
from __future__ import absolute_import

from aiida import orm
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation
from aiida_quantumespresso.data.forceconstants import ForceconstantsData


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

    _internal_retrieve_list = [
        _PHONON_FREQUENCIES_NAME,
        _PHONON_DOS_NAME
    ]
    _default_parser = 'quantumespresso.matdyn'

    @classmethod
    def define(cls, spec):
        super(MatdynCalculation, cls).define(spec)
        spec.input('parent_folder', valid_type=ForceconstantsData, required=True)
        spec.input('kpoints', valid_type=orm.KpointsData, help='Kpoints on which to calculate the phonon frequencies.')
        spec.output('output_parameters', valid_type=orm.Dict)
        spec.output('output_phonon_bands', valid_type=orm.BandsData)
        spec.default_output_node = 'output_parameters'
        spec.exit_code(100, 'ERROR_NO_RETRIEVED_FOLDER',
            message='The retrieved folder data node could not be accessed.')
        spec.exit_code(110, 'ERROR_READING_OUTPUT_FILE',
            message='The output file could not be read from the retrieved folder.')
        spec.exit_code(130, 'ERROR_JOB_NOT_DONE',
            message='The computation did not finish properly ("JOB DONE" not found).')
        spec.exit_code(131, 'ERROR_OUTPUT_KPOINTS_MISSING',
            message='Number of kpoints not found in the output data')
        spec.exit_code(132, 'ERROR_OUTPUT_KPOINTS_INCOMMENSURATE',
            message='Number of kpoints in the inputs is not commensurate with those in the output')

    def _get_following_text(self):
        """Add the kpoints after the namelist."""
        try:
            kpoints = self.inputs.kpoints.get_kpoints()
        except AttributeError:
            kpoints = self.inputs.kpoints.get_kpoints_mesh(print_list=True)

        kpoints_string = [u'{}'.format(len(kpoints))]
        for kpoint in kpoints:
            kpoints_string.append(u'{:18.10f} {:18.10f} {:18.10f}'.format(*kpoint))

        return '\n'.join(kpoints_string) + '\n'

    def prepare_for_submission(self, folder):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        parent_folder = self.inputs.parent_folder

        self._blocked_keywords.append(('INPUT', 'flfrc', parent_folder.filename))

        return super(MatdynCalculation, self).prepare_for_submission(folder)
