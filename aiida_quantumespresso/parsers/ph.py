# -*- coding: utf-8 -*-
from aiida.parsers.parser import Parser
from aiida.orm.nodes.data.folder import FolderData
from aiida.orm.nodes.data.dict import Dict
from aiida_quantumespresso.parsers import QEOutputParsingError
from aiida_quantumespresso.parsers.raw_parser_ph import parse_raw_ph_output
from aiida_quantumespresso.calculations.ph import PhCalculation

class PhParser(Parser):
    """
    This class is the implementation of the Parser class for PHonon.
    """

    def __init__(self,calculation):
        """
        Initialize the instance of PhParser
        """
        # check for valid input
        if not isinstance(calculation,PhCalculation):
            raise QEOutputParsingError("Input calc must be a PhCalculation")
        
        self._calc = calculation
        
        super(PhParser, self).__init__(calculation)
        
    def parse_with_retrieved(self,retrieved):
        """
        Receives in input a dictionary of retrieved nodes.
        Does all the logic here.
        """       
        from aiida.common.exceptions import InvalidOperation
        import os
        
        successful = True
        
        # Check that the retrieved folder is there 
        try:
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
        except KeyError:
            self.logger.error("No retrieved folder found")
            return False, ()

        # check what is inside the folder
        list_of_files = out_folder.get_folder_list()
        # at least the stdout should exist
        if not self._calc._OUTPUT_FILE_NAME in list_of_files:
            successful = False
            new_nodes_tuple = ()
            self.logger.error("Standard output not found")
            return successful,new_nodes_tuple
        
        # if there is something more, I note it down, so to call the raw parser
        # with the right options
        # look for xml
        xml_tensor_file = None
        if self._calc._OUTPUT_XML_TENSOR_FILE_NAME in list_of_files:
            xml_tensor_file = out_folder.get_abs_path(
                                    self._calc._OUTPUT_XML_TENSOR_FILE_NAME)

        # look for dynamical matrices
        dynmat_dir = out_folder.get_abs_path(
                               self._calc._FOLDER_DYNAMICAL_MATRIX)
        dynamical_matrix_list = [ os.path.join(dynmat_dir,fil) 
                                 for fil in os.listdir(dynmat_dir) 
                                 if os.path.isfile(os.path.join(dynmat_dir,fil) ) 
                                 and fil.startswith( os.path.split(self._calc._OUTPUT_DYNAMICAL_MATRIX_PREFIX)[1] )
                                 and not fil.endswith(".freq") ]
        
        # sort according to the number at the end of the dynamical matrix files,
        # when there is one
        try:
            dynamical_matrix_list = sorted( dynamical_matrix_list,
                                        key=lambda x: int(x.split('-')[-1]) )
        except ValueError:
            dynamical_matrix_list = sorted( dynamical_matrix_list)
        
        # define output file name
        out_file = out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME)
        
        # call the raw parsing function
        out_dict,raw_successful = parse_raw_ph_output(out_file,xml_tensor_file,
                                                      dynamical_matrix_list)
        successful = raw_successful if successful else successful
        
        # convert the dictionary into an AiiDA object
        output_params = Dict(dict=out_dict)
        
        # save it into db
        new_nodes_list = [ (self.get_linkname_outparams(),output_params) ]
        
        return successful,new_nodes_list
        
