import numpy as np
from aiida import orm
from qe_tools import CONSTANTS

from aiida_quantumespresso.calculations import _uppercase_dict
from aiida_quantumespresso.calculations.matdyn import MatdynCalculation
from aiida_quantumespresso.utils.mapping import get_logging_container

from .base import BaseParser


class MatdynParser(BaseParser):
    """``Parser`` implementation for the ``MatDynCalculation`` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files from a ``MatdynCalculation`` into output nodes."""
        logs = get_logging_container()

        _, parsed_data, logs = self.parse_stdout_from_retrieved(logs)

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        self.out('output_parameters', orm.Dict(parsed_data))

        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE' in logs.error:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE, logs)

        filename_frequencies = MatdynCalculation._PHONON_FREQUENCIES_NAME  # noqa: SLF001
        filename_dos = MatdynCalculation._PHONON_DOS_NAME  # noqa: SLF001

        if filename_frequencies not in self.retrieved.base.repository.list_object_names():
            return self.exit(self.exit_codes.ERROR_OUTPUT_FREQUENCIES)

        # Extract the kpoints from the input data and create the `KpointsData` for the `BandsData`
        try:
            kpoints = self.node.inputs.kpoints.get_kpoints()
            kpoints_for_bands = self.node.inputs.kpoints.clone()
        except AttributeError:
            kpoints = self.node.inputs.kpoints.get_kpoints_mesh(print_list=True)
            kpoints_for_bands = orm.KpointsData()
            kpoints_for_bands.set_kpoints(kpoints)

        parsed_data = parse_raw_matdyn_phonon_file(
            self.retrieved.base.repository.get_object_content(filename_frequencies)
        )

        if 'parameters' in self.node.inputs:
            parameters = _uppercase_dict(self.node.inputs.parameters.get_dict(), dict_name='parameters')
        else:
            parameters = {'INPUT': {}}

        if parameters['INPUT'].get('dos', False):
            if filename_dos not in self.retrieved.base.repository.list_object_names():
                return self.exit(self.exit_codes.ERROR_OUTPUT_DOS)

            parsed_data.pop('phonon_bands', None)

            with self.retrieved.open(filename_dos) as handle:
                dos_array = np.genfromtxt(handle)

            output_dos = orm.XyData()
            output_dos.set_x(dos_array[:, 0], 'frequency', 'cm^(-1)')
            output_dos.set_y(dos_array[:, 1], 'dos', 'states * cm')

            self.out('output_phonon_dos', output_dos)

        else:
            try:
                num_kpoints = parsed_data.pop('num_kpoints')
            except KeyError:
                return self.exit(self.exit_codes.ERROR_OUTPUT_KPOINTS_MISSING)

            if num_kpoints != kpoints.shape[0]:
                return self.exit(self.exit_codes.ERROR_OUTPUT_KPOINTS_INCOMMENSURATE)

            output_bands = orm.BandsData()
            output_bands.set_kpointsdata(kpoints_for_bands)
            output_bands.set_bands(parsed_data.pop('phonon_bands'), units='THz')

            self.out('output_phonon_bands', output_bands)

        for message in parsed_data['warnings']:
            self.logger.error(message)

        return self.exit(logs=logs)


def parse_raw_matdyn_phonon_file(phonon_frequencies: str) -> dict:
    """Parses the phonon frequencies file.

    :param phonon_frequencies: phonon frequencies file from the matdyn calculation

    :return dict parsed_data: keys:
         * warnings: parser warnings raised
         * num_kpoints: number of kpoints read from the file
         * phonon_bands: BandsData object with the bands for each kpoint
    """
    import re
    import numpy as np

    parsed_data = {}
    parsed_data['warnings'] = []

    lines = phonon_frequencies.splitlines()

    # extract number of bands and kpoints from the header
    # example header line: " &plot nbnd=   6, nks=   1 /"
    header_pattern = re.compile(r'\s*&plot\s+nbnd=\s*(\d+),\s+nks=\s*(\d+)\s*/')
    header_match = re.match(header_pattern, lines.pop(0))
    if not header_match:
        parsed_data['warnings'].append('Number of bands or kpoints unreadable in phonon frequencies file')
        return parsed_data
    num_bands = int(header_match.group(1))
    num_kpoints = int(header_match.group(2))
    parsed_data['num_kpoints'] = num_kpoints

    # initialize array of frequencies
    freq_matrix = np.zeros((num_kpoints, num_bands))

    # In the file, each kpoint block consists of:
    # 1 line with kpoint coordinates (optionally followed by weight)
    # one or more lines with frequencies
    # (maybe up to 6 frequencies per line but it can vary so won't assume that)

    # The blocks will be processed in a loop over the number of kpoints
    # and frequencies will be gradually extracted until the expected number of bands is reached.

    # regex patterns
    # kpoint line ex: "            0.000000  0.000000  0.000000"
    # or with weight: "            0.500000  0.288675  0.000000  0.000000"
    kpoint_pattern = re.compile(r'\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)(?:\s+([-+]?\d+\.\d+))?')
    # frequency line ex: " -148.6347 -124.2795   46.3694  100.8722  110.9098  132.1670"
    # or with attached signs: " -148.70828-124.2696   46.2846  100.8707  110.9253  132.1867"
    frequency_pattern = re.compile(r'\s*([-+]?\d+\.\d+)')

    for kpt_index in range(num_kpoints):
        if not lines:
            parsed_data['warnings'].append('Unexpected end of file while reading kpoints')
            return parsed_data

        kpt_line = lines.pop(0)
        if not re.match(kpoint_pattern, kpt_line):
            parsed_data['warnings'].append(f'Invalid kpoint line format: "{kpt_line}"')
            return parsed_data

        freq_count = 0
        while freq_count < num_bands:
            if not lines:
                parsed_data['warnings'].append('Unexpected end of file while reading frequencies')
                return parsed_data

            freq_line = lines.pop(0)
            freq_matches = re.findall(frequency_pattern, freq_line)
            if not freq_matches:
                parsed_data['warnings'].append(f'Invalid frequency line format: "{freq_line}"')
                return parsed_data

            for freq_str in freq_matches:
                if freq_count < num_bands:
                    try:
                        freq_matrix[kpt_index, freq_count] = (
                            float(freq_str) * CONSTANTS.invcm_to_THz
                        )  # from cm-1 to THz
                    except ValueError:
                        parsed_data['warnings'].append('Error while parsing the frequencies')
                        return parsed_data
                    freq_count += 1
                else:
                    parsed_data['warnings'].append('More frequencies than expected for a kpoint')
                    break

    parsed_data['phonon_bands'] = freq_matrix

    return parsed_data
