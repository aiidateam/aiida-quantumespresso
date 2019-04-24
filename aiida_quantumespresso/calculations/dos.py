# -*- coding: utf-8 -*-
from __future__ import absolute_import
from aiida.common import InputValidationError
from aiida.orm import RemoteData, FolderData, CalcJobNode
from aiida.orm import Dict, XyData
from aiida_quantumespresso.calculations.namelists import NamelistsCalculation
from aiida_quantumespresso.calculations.pw import PwCalculation
import six


class DosCalculation(NamelistsCalculation):
    """
    Plugin for the dos.x code of the Quantum ESPRESSO distribution. Handles
    density of states calculations, and stores the resulting dos arrays and
    integrated dos arrays.
    For more information regarding dos.x
    refer to http://www.quantum-espresso.org/
    """
    
    _DOS_FILENAME = 'aiida.dos'
    _default_namelists = ['DOS']
    _blocked_keywords = [
        ('DOS', 'fildos', _DOS_FILENAME),
        ('DOS', 'outdir', NamelistsCalculation._OUTPUT_SUBFOLDER),
        ('DOS', 'prefix', NamelistsCalculation._PREFIX),
    ]
    _internal_retrieve_list = [_DOS_FILENAME]

    @classmethod
    def define(cls, spec):
        super(DosCalculation, cls).define(spec)
        spec.input('metadata.options.input_filename', valid_type=six.string_types, default=cls._DEFAULT_INPUT_FILE, non_db=True)
        spec.input('metadata.options.output_filename', valid_type=six.string_types, default=cls._DEFAULT_OUTPUT_FILE, non_db=True)
        spec.input('metadata.options.parser_name', valid_type=six.string_types, default='quantumespresso.dos', non_db=True)
        spec.input('parent_folder', valid_type=(RemoteData, FolderData), required=False)
        spec.output('output_parameters', valid_type=Dict)
        spec.output('output_dos', valid_type=XyData)
        spec.exit_code(
            100, 'ERROR_NO_RETRIEVED_FOLDER', message='The retrieved folder data node could not be accessed.')
        spec.exit_code(
            110, 'ERROR_READING_OUTPUT_FILE', message='The output file could not be read from the retrieved folder.')
        spec.exit_code(
            111, 'ERROR_READING_DOS_FILE', message='The dos file could not be read from the retrieved folder.')
        spec.exit_code(
            130, 'ERROR_JOB_NOT_DONE', message='The computation did not finish properly (\'JOB DONE\' not found).')
    
    # TODO: if you want to be able to specify an input calculation (rather than a folder),
    #       you must add this as an input to the NamelistsCalculation with non_db=True, like this:
    # spec.input('parent_calculation', valid_type=CalcJobNode, required=False, non_db=True)
    # TODO: Then you must add the following logic to the prepare_for_submission method
    #       (note that you can't actually modify self.inputs, because it's a 'AttributesFrozendict')
    # def prepare_for_submission(self, folder):
    #     """
    #     Set the parent folder if necessary, or raise if neither parent_folder nor parent_calculation were given.
    #     Then call define() from the super class.
    #     """
    #     if 'parent_folder' in self.inputs and 'parent_calculation' in self.inputs:
    #         raise InputValidationError('You cannot specify both inputs \'parent_folder\' and \'parent_calculation\'.')
    #     elif 'parent_folder' in self.inputs:
    #         pass
    #     elif 'parent_calculation' in self.inputs:
    #         parent_calc = self.inputs.parent_calculation
    #         try:
    #             # TODO: check that the calc's process_type is quantumespresso.pw 
    #             parent_folder = parent_calc.get_outgoing(node_class=RemoteData, link_label_filter='remote_folder').one().node
    #         except ValueError:
    #             raise InputValidationError('Parent calculation does not have a remote folder output node.')
    #         self.inputs['parent_folder'] = parent_folder
    #     else:
    #         raise InputValidationError('Either \'parent_calculation\' or \'parent_folder\' must be in the inputs.')
    # 
    #     super(DosCalculation, self).prepare_for_submission(folder)
    # 
    # TODO: ... except that it won't work, because I can't validate the process_type of the parent calc from the
    #  NamelistsCalculation class: the acceptable types may change according to what is the concrete
    #  type (e.g. DosCalculation, ...).    