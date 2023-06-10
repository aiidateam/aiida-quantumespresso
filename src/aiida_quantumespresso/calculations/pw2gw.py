# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the pw2gw.x code of Quantum ESPRESSO."""
import pathlib

from aiida import orm
from aiida.common import datastructures

from aiida_quantumespresso.calculations.namelists import NamelistsCalculation


class Pw2gwCalculation(NamelistsCalculation):
    """`CalcJob` implementation for the pw2gw.x code of Quantum ESPRESSO."""

    _INPUT_PSEUDOFOLDER = './pseudo/'

    _EPS_X = 'epsX.dat'
    _EPS_Y = 'epsY.dat'
    _EPS_Z = 'epsZ.dat'
    _EPS_TOT = 'epsTOT.dat'

    _default_namelists = ['INPUTPP']
    _internal_retrieve_list = [_EPS_X, _EPS_Y, _EPS_Z, _EPS_TOT]
    _blocked_keywords = [
        ('INPUTPP', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),
        ('INPUTPP', 'prefix', NamelistsCalculation._PREFIX),
    ]

    @classmethod
    def define(cls, spec):
        """Define the process specification."""
        # yapf: disable
        super().define(spec)
        spec.input('parent_folder', valid_type=orm.RemoteData,
            help='Output folder of a completed `PwCalculation`')

        spec.output('output_parameters', valid_type=orm.Dict,
            help='The `output_parameters` output node of the successful calculation.`')
        spec.output('eps', valid_type=orm.ArrayData,
            help='The `eps` output node containing 5 arrays `energy`, `epsX`, `epsY`, `epsZ`, `epsTOT`')

        spec.exit_code(305, 'ERROR_OUTPUT_FILES',
            message='The eps*.dat output files could not be read or parsed.')
        spec.exit_code(330, 'ERROR_OUTPUT_FILES_INVALID_FORMAT',
            message='The eps*.dat output files do not have the expected shape (N, 2).')
        spec.exit_code(331, 'ERROR_OUTPUT_FILES_ENERGY_MISMATCH',
            message='The eps*.dat output files contains different values of energies.')
        spec.exit_code(350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION',
            message='The parser raised an unexpected exception: {exception}')
        # yapf: enable

    def prepare_for_submission(self, folder):
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :class:`~aiida.common.datastructures.CalcInfo` instance.
        """
        calcinfo = super().prepare_for_submission(folder)

        calcinfo.codes_run_mode = datastructures.CodeRunMode.SERIAL

        parent_calc_folder = self.inputs.parent_folder
        calcinfo.remote_copy_list.append((
            parent_calc_folder.computer.uuid,
            str(pathlib.Path(parent_calc_folder.get_remote_path()) / self._INPUT_PSEUDOFOLDER), self._INPUT_PSEUDOFOLDER
        ))

        return calcinfo
