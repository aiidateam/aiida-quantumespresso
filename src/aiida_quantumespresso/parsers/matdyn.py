# -*- coding: utf-8 -*-
from aiida import orm
from qe_tools import CONSTANTS

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

        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE'in logs.error:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE, logs)

        filename_frequencies = MatdynCalculation._PHONON_FREQUENCIES_NAME

        # Extract the kpoints from the input data and create the `KpointsData` for the `BandsData`
        try:
            kpoints = self.node.inputs.kpoints.get_kpoints()
            kpoints_for_bands = self.node.inputs.kpoints.clone()
        except AttributeError:
            kpoints = self.node.inputs.kpoints.get_kpoints_mesh(print_list=True)
            kpoints_for_bands = orm.KpointsData()
            kpoints_for_bands.set_kpoints(kpoints)

        parsed_data = parse_raw_matdyn_phonon_file(self.retrieved.base.repository.get_object_content(filename_frequencies))

        try:
            num_kpoints = parsed_data.pop('num_kpoints')
        except KeyError:
            return self.exit(self.exit_codes.ERROR_OUTPUT_KPOINTS_MISSING, logs)

        if num_kpoints != kpoints.shape[0]:
            return self.exit(self.exit_codes.ERROR_OUTPUT_KPOINTS_INCOMMENSURATE, logs)

        output_bands = orm.BandsData()
        output_bands.set_kpointsdata(kpoints_for_bands)
        output_bands.set_bands(parsed_data.pop('phonon_bands'), units='THz')

        for message in parsed_data['warnings']:
            self.logger.error(message)

        self.out('output_phonon_bands', output_bands)

        return self.exit(logs=logs)


def parse_raw_matdyn_phonon_file(phonon_frequencies):
    """Parses the phonon frequencies file.

    :param phonon_frequencies: phonon frequencies file from the matdyn calculation

    :return dict parsed_data: keys:
         * warnings: parser warnings raised
         * num_kpoints: number of kpoints read from the file
         * phonon_bands: BandsData object with the bands for each kpoint
    """
    import re

    import numpy

    parsed_data = {}
    parsed_data['warnings'] = []

    # extract numbere of bands and kpoints
    try:
        num_bands = int(phonon_frequencies.split('=')[1].split(',')[0])
        num_kpoints = int(phonon_frequencies.split('=')[2].split('/')[0])
        parsed_data['num_kpoints'] = num_kpoints
    except (ValueError, IndexError):
        parsed_data['warnings'].append('Number of bands or kpoints unreadable in phonon frequencies file')
        return parsed_data

    # initialize array of frequencies
    freq_matrix = numpy.zeros((num_kpoints, num_bands))

    split_data = phonon_frequencies.split()
    # discard the header of the file
    raw_data = split_data[split_data.index('/') + 1:]

    # try to improve matdyn deficiencies
    corrected_data = []
    for b in raw_data:
        try:
            corrected_data.append(float(b))
        except ValueError:
            # case in which there are two frequencies attached like -1204.1234-1020.536
            if '-' in b:
                c = re.split('(-)', b)
                d = [i for i in c if i != '']
                for i in range(0, len(d), 2):  # d should have an even number of elements
                    corrected_data.append(float(d[i] + d[i + 1]))
            else:
                # I don't know what to do
                parsed_data['warnings'].append('Bad formatting of frequencies')
                return parsed_data

    counter = 3
    for i in range(num_kpoints):
        for j in range(num_bands):
            try:
                freq_matrix[i, j] = corrected_data[counter] * CONSTANTS.invcm_to_THz  # from cm-1 to THz
            except ValueError:
                parsed_data['warnings'].append('Error while parsing the frequencies')
            except IndexError:
                parsed_data['warnings'].append('Error while parsing the frequencies, dimension exceeded')
                return parsed_data
            counter += 1
        counter += 3  # move past the kpoint coordinates

    parsed_data['phonon_bands'] = freq_matrix

    return parsed_data
