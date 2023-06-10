# -*- coding: utf-8 -*-
from typing import Tuple

from aiida.common import AttributeDict
from aiida.orm import Dict, KpointsData

from aiida_quantumespresso.parsers.base import BaseParser
from aiida_quantumespresso.utils.mapping import get_logging_container


class OpenGridParser(BaseParser):
    """``Parser`` implementation for the ``OpenGridCalculation`` calculation job class."""

    class_error_map = {
        'incompatible FFT grid': 'ERROR_INCOMPATIBLE_FFT_GRID'
    }

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed ``OpenGridCalculation`` into output nodes."""
        logs = get_logging_container()

        stdout, parsed_data, logs = self.parse_stdout_from_retrieved(logs)

        base_exit_code = self.check_base_errors(logs)
        if base_exit_code:
            return self.exit(base_exit_code, logs)

        self.out('output_parameters', Dict(parsed_data))

        for exit_code in self.get_error_map().values():
            if exit_code in logs.error:
                return self.exit(self.exit_codes.get(exit_code), logs)

        kpoints_mesh, kpoints, logs = self.parse_kpoints(stdout, logs)

        self.out('kpoints_mesh', kpoints_mesh)
        self.out('kpoints', kpoints)

        for exit_code in ['ERROR_OUTPUT_KPOINTS_MISMATCH', 'ERROR_OUTPUT_STDOUT_INCOMPLETE']:
            if exit_code in logs.error:
                return self.exit(self.exit_codes.get(exit_code), logs)

        return self.exit(logs=logs)

    @staticmethod
    def parse_kpoints(stdout: str, logs: AttributeDict) -> Tuple[KpointsData, KpointsData, AttributeDict]:
        """Parse the ``stdout`` for the mesh and explicit list of kpoints."""
        lines = stdout.split('\n')

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
            logs.error.append('ERROR_OUTPUT_KPOINTS_MISMATCH')

        return kpoints_mesh, kpoints_list, logs
