# -*- coding: utf-8 -*-
"""Classes and methods for running the ufp2plotcore.sh shell script using AiiDA."""
import os
from aiida import orm
from aiida.common import exceptions
from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida.common.folders import Folder
from aiida.engine import CalcJob, CalcJobProcessSpec
from aiida.plugins import DataFactory

pseudo_upf = DataFactory('pseudo.upf')
upf = DataFactory('upf')
SingleFileData = DataFactory('singlefile')

class UPF2PlotcoreCalculation(CalcJob):
    """Calcjob class to run the ufp2plotcore shell script."""
    
    _default_parser = 'quantumespresso.upf2plotcore'

    @classmethod
    def define(cls, spec: CalcJobProcessSpec):
        """Define the process specification."""
        
        super().define(spec)
        spec.input('upf_data', valid_type=(pseudo_upf, upf), required=True, help='UPF data node to be read by upf2plotcore.sh. This must be an aiida-pseudo type of node ("pseudo.upf") and contain GIPAW information in order to return any meaningful data')
        spec.output('core_wfc_data', valid_type=SingleFileData)
        
        spec.inputs['metadata']['options']['input_filename'].default = 'plotcore.in'
        spec.inputs['metadata']['options']['output_filename'].default = 'plotcore.out'
        spec.inputs['metadata']['options']['parser_name'].default = 'quantumespresso.upf2plotcore'
        spec.inputs['metadata']['options']['resources'].default = {'num_machines': 1, 'num_mpiprocs_per_machine': 1}
        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ',
            message='The stdout output file could not be read.')
        spec.exit_code(315, 'ERROR_OUTPUT_EMPTY', message='No data were retrieved from the pseudopotential. Please ensure that the chosen upf file contains GIPAW information')
        
    def prepare_for_submission(self, folder: Folder) -> CalcInfo:
        """Prepare the calculation for submission. 
        
        Convert the input nodes into the corresponding input files in the format that the 
        code will expect. In addition, define and return a `CalcInfo` instance, which is a 
        simple data structure that contains information for the engine, for example, on what 
        files to copy to the remote machine, what files to retrieve once it has completed,
        specific scheduler settings and more.
        
        :param folder: a temporary folder on the local file system
        :returns: the CalcInfo instance
        """
        
        # upf2plotcore.sh expects a UPF-format (v1 or v2) plaintext file as input
        with folder.open(self.options.input_filename, 'w', encoding='utf8') as handle:
            handle.write(self.inputs.upf_data.get_content())
                
        codeinfo = CodeInfo()
        codeinfo.code_uuid = self.inputs.code.uuid
        codeinfo.stdin_name = self.options.input_filename
        codeinfo.stdout_name = self.options.output_filename
        
        # Create a calcinfo object which contains the information on all codes to use and
        # which outputs to retrieve.
        
        calcinfo = CalcInfo()
        calcinfo.codes_info = [codeinfo]
        calcinfo.retrieve_list = [self.options.output_filename]
        
        return calcinfo
    
