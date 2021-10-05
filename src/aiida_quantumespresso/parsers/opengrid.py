# -*- coding: utf-8 -*-
from aiida.common import NotExistent
from aiida.orm import KpointsData

from aiida_quantumespresso.parsers.parse_raw.base import parse_output_base
from aiida_quantumespresso.parsers.base import Parser


class OpengridParser(Parser):
    """`Parser` implementation for the `OpengridCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed `OpengridCalculation` into output nodes."""
        try:
            out_folder = self.retrieved
        except NotExistent:
            return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_FOLDER)

        try:
            filename_stdout = self.node.get_option('output_filename')
            with out_folder.open(filename_stdout, 'r') as f:
                out_file = f.read()
        except OSError:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_READ)

        parsed_data, logs = parse_output_base(out_file, codename='OPEN_GRID')
        self.emit_logs(logs)

        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE' in logs.error:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE)
        elif logs.error:
            return self.exit(self.exit_codes.ERROR_GENERIC_QE_ERROR)

        try:
            kpoints_mesh, kpoints = self.parse_kpoints(out_file)
        except ValueError:
            return self.exit(self.exit_codes.ERROR_OUTPUT_KPOINTS_MISMATCH)

        # Output both the dimensions and the explict list of kpoints
        self.out('kpoints_mesh', kpoints_mesh)
        self.out('kpoints', kpoints)

    @classmethod
    def parse_kpoints(cls, out_file):
        """Parse and output the dimensions and the explicit list of kpoints"""
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
