import re
from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.base import Parser

class UPF2PlotcoreParser(Parser):
    """ Parser for the UPF2Plotcore calcjob plugin """

    def parse(self, **kwargs):
        """Parse the contents of the output files stored in the `retrieved` output node."""
        from aiida.plugins import DataFactory
        SingleFileData = DataFactory('singlefile')
        
        retrieved = self.retrieved
        try:
            out_file = self.node.get_option('output_filename')
            with retrieved.open(out_file, 'r') as fil:
                lines = fil.readlines()
        except OSError:
            return self.exit(self.exit_codes.ERROR_OUTPUT_STDOUT_READ)
        
        out_file = self.node.get_option('output_filename')
        with retrieved.open(out_file, 'r') as fil:
            lines = fil.readlines()
            if len(lines) == 1:
                return self.exit(self.exit_codes.ERROR_OUTPUT_EMPTY)

        OutFilePath = self.node.get_remote_workdir()
        OutFileName = self.node.get_option('output_filename')
        
        core_wfc_data = OutFilePath + '/' + OutFileName
        
        self.out('core_wfc_data', SingleFileData(core_wfc_data))

