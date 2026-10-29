# -*- coding: utf-8 -*-
"""`CalcJob` implementation for the average.x code of Quantum ESPRESSO."""

from os import path
from aiida import orm
from aiida.common import datastructures
from aiida.engine import CalcJob
from ase.units import Bohr


class AverageCalculation(CalcJob):
    """AiiDA calculation plugin for the `average.x` code of Quantum ESPRESSO.

    This class defines the interface for running the `average.x` post-processing code within the AiiDA framework.
    It specifies the required and optional inputs, expected outputs, and handles the preparation of input files
    and retrieval of output files for the calculation of planar p(z) and macroscopic m(z) averages.

    Inputs:
        - parent_folder (RemoteData or FolderData, required): Output folder of a completed `PpCalculation` where the 'aiida.filplot' file will be copied from.
        - average_axis (Int, optional, default=3): Planar average is done in a plane orthogonal to this axis (1, 2, or 3 for x, y, and z).
        - npts (Int, optional, default=1000): Number of interpolation points for the planar and macroscopic averages along the specified axis.
        - window_size (Float, required): Window size for the macroscopic average in Angstrom (converted to Bohr internally).
        - metadata.options.input_filename (str, default="avg.in"): Name of the input file.
        - metadata.options.output_filename (str, default="avg.out"): Name of the output file.
        - metadata.options.parser_name (str, default="quantumespresso.average"): Name of the parser to use.
        - metadata.options.withmpi (bool, default=False): Whether to run with MPI.

    Outputs:
        - output_parameters (Dict): Parsed parameters from the calculation.
        - output_data (ArrayData): Output data with columns: z in Angstrom, p(z) in eV, m(z) in eV. (Code outputs are in Bohr and Rydberg, conversion is done internally.)

    Exit Codes:
        - 301: ERROR_NO_RETRIEVED_TEMPORARY_FOLDER - The retrieved temporary folder could not be accessed.
        - 302: ERROR_OUTPUT_STDOUT_MISSING - The retrieved folder did not contain the required stdout output file.
        - 303: ERROR_OUTPUT_XML_MISSING - The parent folder did not contain the required XML output file.
        - 310: ERROR_OUTPUT_STDOUT_READ - The stdout output file could not be read.
        - 311: ERROR_OUTPUT_STDOUT_PARSE - The stdout output file could not be parsed.
        - 312: ERROR_OUTPUT_STDOUT_INCOMPLETE - The stdout output file was incomplete.
        - 340: ERROR_OUT_OF_WALLTIME_INTERRUPTED - The calculation stopped prematurely due to walltime exhaustion.
        - 350: ERROR_UNEXPECTED_PARSER_EXCEPTION - The parser raised an unexpected exception.
        - 330: ERROR_OUTPUT_DATAFILE_MISSING - The formatted data output file was not present in the retrieved folder.
        - 331: ERROR_OUTPUT_DATAFILE_READ - The formatted data output file could not be read.
        - 332: ERROR_UNSUPPORTED_DATAFILE_FORMAT - The data file format is not supported by the parser.
        - 333: ERROR_OUTPUT_DATAFILE_PARSE - The formatted data output file could not be parsed.
    """

    _DEFAULT_INPUT_FILE = "avg.in"
    _DEFAULT_OUTPUT_FILE = "avg.out"
    _DEFAULT_INPUT_DATA_FILE = "aiida.filplot"
    _DEFAULT_OUTPUT_DATA_FILE = "avg.dat"

    @classmethod
    def define(cls, spec):
        """Define inputs and outputs of the calculation."""
        super(AverageCalculation, cls).define(spec)

        # Input ports
        spec.input(
            "parent_folder",
            valid_type=(orm.RemoteData, orm.FolderData),
            required=True,
            help="Output folder of a completed `PpCalculation` where the 'aiida.filplot' file will be copied from.",
        )
        spec.input(
            "average_axis",
            valid_type=orm.Int,
            required=False,
            default=lambda: orm.Int(3),
            help="Planar average done in plane orthogonal to this axis (1,2 or 3 for x,y and z)",
        )
        spec.input(
            "npts",
            valid_type=orm.Int,
            required=False,
            default=lambda: orm.Int(1000),
            help="Number of interpolation points of the planar and macroscopic averages.\
                If less points than the nb of FFT points in average_axis direction, no interpolation is done.",
        )
        spec.input(
            "window_size",
            valid_type=orm.Float,
            required=True,
            help="Window size for the macroscopic average in Angstrom (the code reads in Bohr, conversion is done internally)",
        )

        spec.input(
            "metadata.options.input_filename",
            valid_type=str,
            default=cls._DEFAULT_INPUT_FILE,
        )
        spec.input(
            "metadata.options.output_filename",
            valid_type=str,
            default=cls._DEFAULT_OUTPUT_FILE,
        )
        spec.input(
            "metadata.options.parser_name",
            valid_type=str,
            default="quantumespresso.average",
        )
        spec.input("metadata.options.withmpi", valid_type=bool, default=False)

        # Output ports
        spec.output(
            "output_parameters",
            valid_type=orm.Dict,
        )
        spec.output(
            "output_data",
            valid_type=orm.ArrayData,
            help="The output data with columns: \n\
z in Ang, p(z) in eV, m(z) in eV. Code outputs are in Bohr and Rydberg. Conversion is done internally.",
        )

        # Standard exceptions
        spec.exit_code(
            301,
            "ERROR_NO_RETRIEVED_TEMPORARY_FOLDER",
            message="The retrieved temporary folder could not be accessed.",
        )
        spec.exit_code(
            302,
            "ERROR_OUTPUT_STDOUT_MISSING",
            message="The retrieved folder did not contain the required stdout output file.",
        )
        spec.exit_code(
            303,
            "ERROR_OUTPUT_XML_MISSING",
            message="The parent folder did not contain the required XML output file.",
        )
        spec.exit_code(
            310,
            "ERROR_OUTPUT_STDOUT_READ",
            message="The stdout output file could not be read.",
        )
        spec.exit_code(
            311,
            "ERROR_OUTPUT_STDOUT_PARSE",
            message="The stdout output file could not be parsed.",
        )
        spec.exit_code(
            312,
            "ERROR_OUTPUT_STDOUT_INCOMPLETE",
            message="The stdout output file was incomplete.",
        )
        spec.exit_code(
            340,
            "ERROR_OUT_OF_WALLTIME_INTERRUPTED",
            message="The calculation stopped prematurely because it ran out of walltime but the job was killed by the "
            "scheduler before the files were safely written to disk for a potential restart.",
        )
        spec.exit_code(
            350,
            "ERROR_UNEXPECTED_PARSER_EXCEPTION",
            message="The parser raised an unexpected exception: {exception}",
        )

        # Output datafile related exceptions
        spec.exit_code(
            330,
            "ERROR_OUTPUT_DATAFILE_MISSING",
            message="The formatted data output file `{filename}` was not present in the retrieved (temporary) folder.",
        )
        spec.exit_code(
            331,
            "ERROR_OUTPUT_DATAFILE_READ",
            message="The formatted data output file `{filename}` could not be read.",
        )
        spec.exit_code(
            332,
            "ERROR_UNSUPPORTED_DATAFILE_FORMAT",
            message="The data file format is not supported by the parser",
        )
        spec.exit_code(
            333,
            "ERROR_OUTPUT_DATAFILE_PARSE",
            message="The formatted data output file `{filename}` could not be parsed: {exception}",
        )

    def prepare_for_submission(
        self, folder
    ):  # pylint: disable=too-many-branches,too-many-statements
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :class:`~aiida.common.datastructures.CalcInfo` instance.
        """
        parent_folder = self.inputs.parent_folder

        remote_copy_list = [
            parent_folder.computer.uuid,
            path.join(
                parent_folder.get_remote_path(),
                self._DEFAULT_INPUT_DATA_FILE,
            ),
            self._DEFAULT_INPUT_DATA_FILE,
        ]
        local_copy_list = []

        # Code information
        codeinfo = datastructures.CodeInfo()
        codeinfo.stdin_name = self.inputs.metadata.options.input_filename
        codeinfo.stdout_name = self.inputs.metadata.options.output_filename
        codeinfo.code_uuid = self.inputs.code.uuid

        # Calculation information
        calcinfo = datastructures.CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.remote_copy_list = remote_copy_list
        calcinfo.local_copy_list = local_copy_list

        # Retrieve by default the output file and the data file
        calcinfo.retrieve_list = [self.inputs.metadata.options.output_filename]
        calcinfo.retrieve_temporary_list = [self._DEFAULT_OUTPUT_DATA_FILE]

        # Create the input file
        input_filename = self.inputs.metadata.options.input_filename
        with folder.open(input_filename, "w") as infile:
            infile.write("1\n")
            infile.write(self._DEFAULT_INPUT_DATA_FILE + "\n")
            infile.write("1.D0\n")
            infile.write(f"{self.inputs.npts.value}\n")
            infile.write(f"{self.inputs.average_axis.value}\n")
            infile.write(f"{self.inputs.window_size.value / Bohr:.8f}\n")

        return calcinfo
