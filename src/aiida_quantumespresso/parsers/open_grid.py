# -*- coding: utf-8 -*-
from aiida.orm import Dict, KpointsData

from aiida_quantumespresso.parsers.base import BaseParser


class OpenGridParser(BaseParser):
    """``Parser`` implementation for the ``OpenGridCalculation`` calculation job class."""

    class_error_map = {
        'incompatible FFT grid': 'ERROR_INCOMPATIBLE_FFT_GRID'
    }

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed ``OpenGridCalculation`` into output nodes."""
        parsed_stdout, logs_stdout = self._parse_stdout_from_retrieved()
        self._emit_logs(logs_stdout)

        kpoints_mesh = parsed_stdout.pop('kpoints_mesh', None)
        kpoints = parsed_stdout.pop('kpoints', None)

        self.out('output_parameters', Dict(parsed_stdout))

        for exit_code in [
            'ERROR_INCOMPATIBLE_FFT_GRID', 'ERROR_OUTPUT_KPOINTS_MISMATCH', 'ERROR_OUTPUT_STDOUT_MISSING',
            'ERROR_OUTPUT_STDOUT_READ', 'ERROR_OUTPUT_STDOUT_INCOMPLETE'
            ]:
            if exit_code in logs_stdout.error:
                return self._exit(self.exit_codes.get(exit_code))

        # Output both the dimensions and the explict list of kpoints
        self.out('kpoints_mesh', kpoints_mesh)
        self.out('kpoints', kpoints)

    @classmethod
    def parse_stdout(cls, stdout: str) -> tuple:
        """Parse the ``stdout`` content of a Quantum ESPRESSO ``open_grid.x`` calculation.

        Overridden from parent method to include kpoints parsing.

        :param stdout: the stdout content as a string.
        :returns: tuple of two dictionaries, with the parsed data and log messages, respectively.
        """
        parsed_data, logs = super().parse_stdout(stdout)

        if logs['error']:
            return parsed_data, logs

        kpoints_mesh, kpoints = cls.parse_kpoints(stdout)

        parsed_data['kpoints_mesh'] = kpoints_mesh
        parsed_data['kpoints'] = kpoints

        return parsed_data, logs

    @staticmethod
    def parse_kpoints(out_file):
        """Parse and output the dimensions and the explicit list of kpoints."""
        lines = out_file.split('\n')

        kpoints = []
        weights = []
        found_kpoints = False

        for line in lines:
            if 'EXX: q-point mesh:' in line:
                kmesh = [int(i) for i in line.strip().split()[-3:]]
            if 'List to be put in the .win file of wannier90' in line:
                found_kpoints = True
                continue
            if found_kpoints:
                line = line.strip()
                if line != '':
                    line = [float(i) for i in line.split()]
                    kpoints.append(line[:-1])
                    weights.append(line[-1])
                else:
                    found_kpoints = False

        kpoints_mesh = KpointsData()
        kpoints_mesh.set_kpoints_mesh(kmesh)
        kpoints_list = KpointsData()
        kpoints_list.set_kpoints(kpoints, cartesian=False, weights=weights)

        if kmesh[0] * kmesh[1] * kmesh[2] != len(kpoints):
            raise ValueError('Mismatch between kmesh dimensions and number of kpoints')

        return kpoints_mesh, kpoints_list
