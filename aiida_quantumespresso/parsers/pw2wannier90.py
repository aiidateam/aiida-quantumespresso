# -*- coding: utf-8 -*-
from aiida.orm.data.parameter import ParameterData
from aiida.parsers.parser import Parser
from aiida_quantumespresso.parsers.raw_parser_simple import parse_qe_simple

class Pw2wannier90Parser(Parser):
    """
    This class is the implementation of the Parser class for pw2wannier90.x
    """
    def parse_with_retrieved(self,retrieved):
        """      
        Parses the datafolder, stores results.
        This parser for this simple code does simply store in the DB a node
        representing the file of forces in real space
        """
        from aiida.common.exceptions import InvalidOperation

        # suppose at the start that the job is successful
        successful = True
        warnings = []

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
            self.logger.error("Standard output not found")
            return successful,()
        
        # check that the file has finished (i.e. JOB DONE is inside the file)
        filpath = out_folder.get_abs_path(self._calc._OUTPUT_FILE_NAME)
        with open(filpath,'r') as fil:
            lines = fil.read()

        successful_raw, out_dict = parse_qe_simple(lines, codename="PW2WANNIER")
        # If any failed, it's failed
        successful = successful and successful_raw

        out_params = ParameterData(dict=out_dict)
        new_nodes_list = [(self.get_linkname_outparams(), out_params)]

        return successful,new_nodes_list

