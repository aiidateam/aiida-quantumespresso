# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the pw2wannier.x code of Quantum ESPRESSO."""
from __future__ import absolute_import
from aiida.orm import RemoteData, FolderData, SinglefileData, Dict

from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class Pw2wannier90Calculation(NamelistsCalculation):
    """`CalcJob` implementation for the pw2wannier.x code of Quantum ESPRESSO.

    For more information, refer to http://www.quantum-espresso.org/ and http://www.wannier.org/
    """
    _default_namelists = ['INPUTPP']
    _SEEDNAME = 'aiida'
    _blocked_keywords = [('INPUTPP', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),
                         ('INPUTPP', 'prefix', NamelistsCalculation._PREFIX), ('INPUTPP', 'seedname', _SEEDNAME)]
    # By default we do not download anything else than aiida.out. One can add the files
    # _SEEDNAME.amn/.nnm/.eig to inputs.settings['ADDITIONAL_RETRIEVE_LIST'] to retrieve them.
    _internal_retrieve_list = []
    _default_parser = 'quantumespresso.pw2wannier90'

    @classmethod
    def define(cls, spec):
        # yapf: disable
        super(Pw2wannier90Calculation, cls).define(spec)
        spec.input('nnkp_file', valid_type=SinglefileData,
                   help='A SinglefileData containing the .nnkp file generated by wannier90.x -pp')
        spec.input('parent_folder', valid_type=(RemoteData, FolderData),
                   help='The output folder of a pw.x calculation')
        spec.output('output_parameters', valid_type=Dict)
        spec.default_output_node = 'output_parameters'
        spec.exit_code(
            100, 'ERROR_NO_RETRIEVED_FOLDER', message='The retrieved folder data node could not be accessed.')
        spec.exit_code(
            110, 'ERROR_READING_OUTPUT_FILE', message='The output file could not be read from the retrieved folder.')
        spec.exit_code(
            130, 'ERROR_JOB_NOT_DONE', message='The computation did not finish properly (\'JOB DONE\' not found).')
        spec.exit_code(
            140, 'ERROR_GENERIC_QE_ERROR', message='QE printed an error message')
        spec.exit_code(
            150, 'ERROR_GENERIC_PARSING_FAILURE', message='An error happened while parsing the output file')

    def prepare_for_submission(self, folder):
        """Prepare the inputs of the calculation and the calcinfo data.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """
        # Run the global namelist logic
        calcinfo = super(Pw2wannier90Calculation, self).prepare_for_submission(folder)

        # Put the nnkp in the folder, with the correct filename
        nnkp_file = self.inputs.nnkp_file
        calcinfo.local_copy_list.append(
            (nnkp_file.uuid, nnkp_file.filename, '{}.nnkp'.format(self._SEEDNAME))
        )

        return calcinfo
