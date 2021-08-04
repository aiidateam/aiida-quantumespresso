# -*- coding: utf-8 -*-
from aiida.common import NotExistent, LinkType
from aiida.orm import Dict, KpointsData, CalcJobNode

from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.parse_raw.base import parse_output_base
from aiida_quantumespresso.parsers.base import Parser


class OpengridParser(Parser):
    """`Parser` implementation for the `OpengridCalculation` calculation job class."""

    def parse(self, **kwargs):
        """Parse the retrieved files of a completed `OpengridCalculation` into output nodes.

        Two nodes that are expected are the default 'retrieved' `FolderData` node which will store the retrieved files
        permanently in the repository.
        """
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

        self.parse_kpoints(out_file)

        try:
            parent_calc = (
                self.node.inputs.parent_folder.get_incoming(node_class=CalcJobNode,
                                                            link_type=LinkType.CREATE).one().node
            )
        except ValueError as e:
            raise QEOutputParsingError(f'Could not get parent calculation of input folder: {e}')

        # We need to output additional nodes for a subsequent projwfc calculation:
        # projwfc parser needs `number_of_spin_components`
        try:
            parent_param = parent_calc.get_outgoing(link_label_filter='output_parameters').one().node
        except ValueError:
            raise QEOutputParsingError('The parent had no output_parameters! Cannot parse from this!')
        nspin = parent_param.get_dict()['number_of_spin_components']
        parsed_data['number_of_spin_components'] = nspin
        self.out('output_parameters', Dict(dict=parsed_data))

    def parse_kpoints(self, out_file):
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
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE)

        self.out('kpoints_mesh', kpoints_mesh)
        self.out('kpoints', kpoints_list)
